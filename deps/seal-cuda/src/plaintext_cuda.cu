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

    // void PlaintextCuda::load(std::istream& stream) {
    //     loadt(stream, &parms_id_);
    //     loadt(stream, &coeff_count_);
    //     loadt(stream, &scale_);
    //     size_t dataSize;
    //     loadt(stream, &dataSize);
    //     util::HostArray<pt_coeff_type> host(dataSize);
    //     stream.read(reinterpret_cast<char*>(host.get()), dataSize * sizeof(pt_coeff_type));
    //     data_.ensure(dataSize);
    //     KernelProvider::copy(data_.get(), host.get(), dataSize);
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
            auto data__ = data_.toHost();
            size_t dataSize = data__.size();
            // savet(stream, &dataSize);
            stream.write(reinterpret_cast<const char *>(&dataSize), sizeof(size_t));
            stream.write(reinterpret_cast<const char*>(data__.begin()), sizeof(pt_coeff_type) * data__.size());
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

    void PlaintextCuda::load(std::istream& stream)
    { 
        using namespace std;
        PlaintextCuda new_data;

        auto old_except_mask = stream.exceptions();
        try
        {
            // Throw exceptions on std::ios_base::badbit and std::ios_base::failbit
            stream.exceptions(std::ios_base::badbit | std::ios_base::failbit);

            ParmsID parms_id{};
            stream.read(reinterpret_cast<char *>(&parms_id), sizeof(ParmsID));

            uint64_t coeff_count64 = 0;
            stream.read(reinterpret_cast<char *>(&coeff_count64), sizeof(uint64_t));

            double scale = 0;
            stream.read(reinterpret_cast<char *>(&scale), sizeof(double));

            // Set the metadata
            new_data.parms_id_ = parms_id;
            new_data.coeff_count_ = static_cast<size_t>(coeff_count64);
            new_data.scale_ = scale;
            
            size_t dataSize;
            stream.read(reinterpret_cast<char *>(&dataSize), sizeof(size_t));
            // new_data.data_.resize(dataSize);
            // stream.read(reinterpret_cast<char*>(new_data.data_.begin()), sizeof(pt_coeff_type) * dataSize);

            util::HostArray<pt_coeff_type> host(dataSize);
            stream.read(reinterpret_cast<char*>(host.get()), dataSize * sizeof(pt_coeff_type));
            data_.ensure(dataSize);
            KernelProvider::copy(data_.get(), host.get(), dataSize);
        }
        catch (const ios_base::failure &)
        {
            stream.exceptions(old_except_mask);
            throw runtime_error("I/O error");
        }
        catch (...)
        {
            stream.exceptions(old_except_mask);
            throw;
        }
        stream.exceptions(old_except_mask);

        swap(*this, new_data);
    }

}