#include "plaintext_cuda.cuh"
#include "serialize.h"

namespace troy {
    

    // void PlaintextCuda::save(std::ostream& stream) const {
    //     savet(stream, &parms_id_);
    //     savet(stream, &coeff_count_);
    //     savet(stream, &scale_);
    //     auto r = data_.toHost();
    //     size_t dataSize = r.size();
    //     savet(stream, &dataSize);
    //     stream.write(reinterpret_cast<char*>(r.begin()), sizeof(pt_coeff_type) * r.size());
    // }

    void PlaintextCuda::save(std::ostream &stream) const
    {
        auto old_except_mask = stream.exceptions();
        try
        {
            // Throw exceptions on std::ios_base::badbit and std::ios_base::failbit
            stream.exceptions(std::ios_base::badbit | std::ios_base::failbit);

            stream.write(reinterpret_cast<const char *>(&parms_id_), sizeof(ParmsID));
            uint64_t coeff_count64 = static_cast<uint64_t>(coeff_count_);
            stream.write(reinterpret_cast<const char *>(&coeff_count64), sizeof(uint64_t));
            stream.write(reinterpret_cast<const char *>(&scale_), sizeof(double));
            // data_.save(stream, compr_mode_type::none);
            auto r = data_.toHost();
            uint64_t dataSize = r.size();
            savet(stream, &dataSize);
            stream.write(reinterpret_cast<char*>(r.begin()), sizeof(pt_coeff_type) * r.size());
        }
        catch (const std::ios_base::failure &)
        {
            stream.exceptions(old_except_mask);
            throw std::runtime_error("I/O error");
        }
        catch (...)
        {
            stream.exceptions(old_except_mask);
            throw;
        }
        stream.exceptions(old_except_mask);
    }

    void PlaintextCuda::load(std::istream& stream) {
        loadt(stream, &parms_id_);
        loadt(stream, &coeff_count_);
        loadt(stream, &scale_);
        size_t dataSize;
        loadt(stream, &dataSize);
        util::HostArray<pt_coeff_type> host(dataSize);
        stream.read(reinterpret_cast<char*>(host.get()), dataSize * sizeof(pt_coeff_type));
        data_.ensure(dataSize);
        KernelProvider::copy(data_.get(), host.get(), dataSize);
    }

}