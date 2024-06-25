#include "../src/troy_cuda.cuh"
#include <vector>
#include <string>
#include <sys/time.h>
#include <cassert>
#include <map>
#include <complex>
#include <iomanip>

using namespace troyn;
using std::vector;
using std::complex;

namespace troytest {
    
    template<typename T> 
    vector<T> vectorAdd(const vector<T> a, const vector<T> b) {
        assert(a.size() == b.size());
        vector<T> ret; ret.reserve(a.size());
        for (size_t i = 0; i < a.size(); i++) ret.push_back(a[i] + b[i]);
        return ret;
    }

    template<typename T> 
    vector<T> vectorMultiply(const vector<T> a, const vector<T> b) {
        assert(a.size() == b.size());
        vector<T> ret; ret.reserve(a.size());
        for (size_t i = 0; i < a.size(); i++) ret.push_back(a[i] * b[i]);
        return ret;
    }

    template<typename T> 
    vector<T> vectorSub(const vector<T> a, const vector<T> b) {
        assert(a.size() == b.size());
        vector<T> ret; ret.reserve(a.size());
        for (size_t i = 0; i < a.size(); i++) ret.push_back(a[i] - b[i]);
        return ret;
    }

    template<typename T> 
    vector<T> vectorNegate(const vector<T> a) {
        vector<T> ret; ret.reserve(a.size());
        for (size_t i = 0; i < a.size(); i++) ret.push_back(-a[i]);
        return ret;
    }

    template<typename T>
    vector<T> vectorRotate(const vector<T>& a, size_t n) {
        vector<T> ret; ret.reserve(a.size());
        ret.insert(ret.end(), a.begin() + n, a.end());
        ret.insert(ret.end(), a.begin(), a.begin() + n);
    }

    template<typename T> 
    inline bool is_equal(const vector<T>& a, const vector<T>& b, double eps = 1e-6) {
        if (a.size() != b.size()) return false;
        for (size_t i = 0; i < a.size(); i++) {
            if (std::abs(a[i] - b[i]) > eps) return false;
        }
        return true;
    }

    template<typename T> 
    inline bool is_zero(const vector<T>& a, double eps = 1e-6) {
        for (size_t i = 0; i < a.size(); i++) {
            if (std::abs(a[i]) > eps) return false;
        }
        return true;
    }

    template<typename T>
    void printVector(const vector<T>& r, bool full = false) {
        std::cout << "[";
        for (size_t i = 0; i < r.size(); i++) {
            if (r.size() > 8 && !full && i == 4) {
                std::cout << " ...";
                i = r.size() - 4;
            }
            if (i!=0) std::cout << ", ";
            std::cout << std::setprecision(3) << std::fixed << (double) r[i];
        }
        std::cout << "]" << std::endl;
    }

    inline std::string pass_str(bool pass) {
        std::stringstream ss;
        // color red if failed, green if passed
        ss << (pass ? "\033[1;32m" : "\033[1;31m") << (pass ? "PASS" : "FAIL") << "\033[0m";
        return ss.str();
    }

    class Timer {
    public:
        std::vector<timeval> times;
        std::vector<double> accumulated; // ms
        std::vector<std::string> names;
        Timer() {}
        long registerTimer(std::string name = "") {
            times.push_back(timeval()); 
            accumulated.push_back(0);
            int ret = times.size() - 1;
            names.push_back(name);
            return ret;
        }
        void tick(long i = 0) {
            if (times.size() < 1) registerTimer();
            assert(i < times.size());
            gettimeofday(&times[i], 0);
        }
        double tock(long i = 0) {
            assert(i < times.size());
            timeval s; gettimeofday(&s, 0);
            auto timeElapsed = (s.tv_sec - times[i].tv_sec) * 1000.0;
            timeElapsed += (s.tv_usec - times[i].tv_usec) / 1000.0;
            accumulated[i] += timeElapsed;
            return accumulated[i];
        }
        
        void clear() {
            times.clear();
            accumulated.clear();
            names.clear();
        }

        std::map<std::string, double> gather(double divisor = 1) {
            std::map<std::string, double> p;
            for (long i=0; i<times.size(); i++) {
                p[names[i]] = accumulated[i] / divisor;
            }
            clear();
            return p;
        }
    };

    class TimeTest {
        
    protected:
        Timer tim;
        Encryptor* encryptor;
        Decryptor* decryptor;
        Evaluator* evaluator;
        SEALContext* context;
        RelinKeys rlk;
        PublicKey pk;
        GaloisKeys gk;
        KeyGenerator* keygen;

    public:
        TimeTest() {
            tim.clear();
            encryptor = nullptr;
            evaluator = nullptr;
            context = nullptr;
            decryptor = nullptr;
        }

        ~TimeTest() {
            if (encryptor) delete encryptor;
            if (evaluator) delete evaluator;
            if (context) delete context;
            if (decryptor) delete decryptor;
            if (keygen) delete keygen;
        }

        virtual Plaintext randomPlaintext() = 0;
        virtual Ciphertext randomCiphertext() = 0;
        // virtual void testEncode() = 0;

        void testEncrypt(int repeatCount = 1000) {
            auto p1 = randomPlaintext();
            Ciphertext c2;
            Plaintext p2;
            auto t1 = tim.registerTimer("Encrypt");
            auto t2 = tim.registerTimer("Decrypt");
            for (int t = 0; t < repeatCount; t++) {
                tim.tick(t1);
                encryptor->encrypt(p1, c2);
                tim.tock(t1);
                tim.tick(t2);
                decryptor->decrypt(c2, p2);
                tim.tock(t2);
            }
            printTimer(tim.gather(repeatCount));
        }

        void printTimer(std::map<std::string, double> r) {
            for (auto& p: r) {
                std::cout << std::setw(25) << std::right << p.first << ":";
                std::cout << std::setw(10) << std::right << std::fixed << std::setprecision(3)
                    << p.second << std::endl;
            }
        }

        void testAdd(int repeatCount = 1000) {
            auto c1 = randomCiphertext();
            auto c2 = randomCiphertext();
            Ciphertext c3;
            auto t1 = tim.registerTimer("Add-assign");
            auto t2 = tim.registerTimer("Add-inplace");
            for (int t = 0; t < repeatCount; t++) {
                tim.tick(t1);
                evaluator->add(c1, c2, c3);
                tim.tock(t1);
                tim.tick(t2);
                evaluator->addInplace(c3, c1);
                tim.tock(t2);
            }
            printTimer(tim.gather(repeatCount));
        }

        void testAddPlain(int repeatCount = 1000) {
            auto c1 = randomCiphertext();
            auto p2 = randomPlaintext();
            Ciphertext c3;
            auto t1 = tim.registerTimer("AddPlain-assign");
            auto t2 = tim.registerTimer("AddPlain-inplace");
            for (int t = 0; t < repeatCount; t++) {
                tim.tick(t1);
                evaluator->addPlain(c1, p2, c3);
                tim.tock(t1);
                tim.tick(t2);
                evaluator->addPlainInplace(c3, p2);
                tim.tock(t2);
            }
            printTimer(tim.gather(repeatCount));
        }

        void testMultiplyPlain(int repeatCount = 1000) {
            auto c1 = randomCiphertext();
            auto p2 = randomPlaintext();
            Ciphertext c3;
            auto t1 = tim.registerTimer("MultiplyPlain-assign");
            auto t2 = tim.registerTimer("MultiplyPlain-inplace");
            for (int t = 0; t < repeatCount; t++) {
                tim.tick(t1);
                evaluator->multiplyPlain(c1, p2, c3);
                tim.tock(t1);
                tim.tick(t2);
                evaluator->multiplyPlainInplace(c3, p2);
                tim.tock(t2);
            }
            printTimer(tim.gather(repeatCount));
        }

        void testSquare(int repeatCount = 1000) {
            auto c1 = randomCiphertext();
            Ciphertext c2;
            Ciphertext c3;
            auto t1 = tim.registerTimer("Square-assign");
            auto t2 = tim.registerTimer("Square-inplace");
            for (int t = 0; t < repeatCount; t++) {
                tim.tick(t1);
                evaluator->square(c1, c2);
                tim.tock(t1);
                c3 = c1;
                tim.tick(t2);
                evaluator->squareInplace(c3);
                tim.tock(t2);
            }
            printTimer(tim.gather(repeatCount));
        }

        void testMemoryPool(int repeatCount = 1000) {
            auto t1 = tim.registerTimer("Preallocate");
            auto t2 = tim.registerTimer("Allocate");
            tim.tick(t1);
            auto c1 = randomCiphertext();
            Ciphertext c2;
            for (int t = 0; t < repeatCount; t++) {
                evaluator->square(c1, c2);
            }
            tim.tock(t1);
            tim.tick(t2);
            for (int t = 0; t < repeatCount; t++) {
                Ciphertext c3;
                evaluator->square(c1, c3);
            }
            tim.tock(t2);
            printTimer(tim.gather(repeatCount));
        }

    };

    class TimeTestCKKS: public TimeTest {

        CKKSEncoder* encoder;
        size_t slotCount;
        int dataBound;
        double delta;
    
    public:

        TimeTestCKKS(size_t polyModulusDegree, vector<int> qs, int dataBound = 1<<6, double delta=static_cast<double>(1<<16)) {
            KernelProvider::initialize();
            slotCount = polyModulusDegree / 2;
            this->dataBound = dataBound;
            this->delta = delta;
            EncryptionParameters parms(SchemeType::ckks);
            parms.setPolyModulusDegree(polyModulusDegree);
            parms.setCoeffModulus(CoeffModulus::Create(polyModulusDegree, qs));
            context = new SEALContext(parms);
            keygen = new KeyGenerator(*context);
            keygen->createPublicKey(pk);
            keygen->createRelinKeys(rlk);
            keygen->createGaloisKeys(gk);
            encoder = new CKKSEncoder(*context);
            encryptor = new Encryptor(*context, pk);
            decryptor = new Decryptor(*context, keygen->secretKey());
            evaluator = new Evaluator(*context);
        }

        ~TimeTestCKKS() {
            if (encoder) delete encoder;
        }
        
        static vector<complex<double>> randomVector(size_t count, int data_bound) {
            vector<complex<double>> input(count, 0.0);
            for (size_t i = 0; i < count; i++)
            {
                input[i] = static_cast<double>(rand() % data_bound);
            }
            return input;
        }

        Plaintext randomPlaintext() override {
            auto p = randomVector(slotCount, dataBound);
            Plaintext ret; encoder->encode(p, delta, ret);
            return std::move(ret);
        }

        Ciphertext randomCiphertext() override {
            auto r = randomPlaintext();
            Ciphertext ret; encryptor->encrypt(r, ret);
            return std::move(ret);
        }

        void testEncode(int repeatCount = 1000) {
            auto m1 = randomVector(slotCount, dataBound);
            auto m2 = randomVector(slotCount, dataBound);
            Plaintext p1;
            auto t1 = tim.registerTimer("Encode");
            auto t2 = tim.registerTimer("Decode");
            for (int t = 0; t < repeatCount; t++) {
                tim.tick(t1);
                encoder->encode(m1, delta, p1);
                tim.tock(t1);
                tim.tick(t2);
                encoder->decode(p1, m2);
                tim.tock(t2);
            }
            printTimer(tim.gather(repeatCount));
        }

        void testMultiplyRescale(int repeatCount = 100) {
            auto c1 = randomCiphertext();
            auto c2 = randomCiphertext();
            Ciphertext c3, c4;
            Ciphertext c5;
            auto t1 = tim.registerTimer("Multiply-assign");
            auto t2 = tim.registerTimer("Rescale-assign");
            auto t3 = tim.registerTimer("Multiply-inplace");
            auto t4 = tim.registerTimer("Rescale-inplace");
            for (int t = 0; t < repeatCount; t++) {
                tim.tick(t1);
                evaluator->multiply(c1, c2, c3);
                tim.tock(t1);
                tim.tick(t2);
                evaluator->rescaleToNext(c3, c4);
                tim.tock(t2);
                c5 = c1;
                tim.tick(t3);
                evaluator->multiplyInplace(c5, c2);
                tim.tock(t3);
                tim.tick(t4);
                evaluator->rescaleToNextInplace(c5);
                tim.tock(t4);
            }
            printTimer(tim.gather(repeatCount));
        }

        void testRotateVector(int repeatCount = 100) {
            auto c1 = randomCiphertext();
            Ciphertext c2;
            auto t1 = tim.registerTimer("Rotate-assign");
            auto t2 = tim.registerTimer("Rotate-inplace");
            for (int t = 0; t < repeatCount; t++) {
                tim.tick(t1);
                evaluator->rotateVector(c1, 1, gk, c2);
                tim.tock(t1);
                tim.tick(t2);
                evaluator->rotateVectorInplace(c1, 1, gk);
                tim.tock(t2);
            }
            printTimer(tim.gather(repeatCount));
        }

        void testAll() {
            this->testEncode();
            this->testEncrypt();
            this->testAdd();
            this->testAddPlain();
            this->testMultiplyRescale();
            this->testMultiplyPlain();
            this->testSquare();
            this->testRotateVector();
            this->testMemoryPool();
        }

    };



    class TimeTestBFVBGV: public TimeTest {

        BatchEncoder* encoder;
        size_t slotCount;
        int dataBound;
        double delta;
    
    public:

        TimeTestBFVBGV(bool bgv, size_t polyModulusDegree, uint64_t plainModulusBitSize, vector<int> qs, int dataBound = 1<<6) {
            KernelProvider::initialize();
            slotCount = polyModulusDegree;
            this->dataBound = dataBound;
            this->delta = delta;
            EncryptionParameters parms(bgv ? SchemeType::bgv : SchemeType::bfv);
            parms.setPolyModulusDegree(polyModulusDegree);
            parms.setPlainModulus(PlainModulus::Batching(polyModulusDegree, plainModulusBitSize));
            // parms.setCoeffModulus(CoeffModulus::BFVDefault(polyModulusDegree));
            parms.setCoeffModulus(CoeffModulus::Create(polyModulusDegree, qs));
            context = new SEALContext(parms);
            keygen = new KeyGenerator(*context);
            keygen->createPublicKey(pk);
            keygen->createRelinKeys(rlk);
            keygen->createGaloisKeys(gk);
            encoder = new BatchEncoder(*context);
            encryptor = new Encryptor(*context, pk);
            decryptor = new Decryptor(*context, keygen->secretKey());
            evaluator = new Evaluator(*context);
        }

        ~TimeTestBFVBGV() {
            if (encoder) delete encoder;
        }
        
        static vector<int64_t> randomVector(size_t count, int data_bound) {
            vector<int64_t> input(count, 0);
            for (size_t i = 0; i < count; i++)
            {
                input[i] = rand() % data_bound;
            }
            return input;
        }

        Plaintext randomPlaintext() override {
            auto p = randomVector(slotCount, dataBound);
            Plaintext ret; encoder->encode(p, ret);
            return std::move(ret);
        }

        Ciphertext randomCiphertext() override {
            auto r = randomPlaintext();
            Ciphertext ret; encryptor->encrypt(r, ret);
            return std::move(ret);
        }

        vector<int64_t> decode(const Plaintext& p) {
            vector<int64_t> r;
            encoder->decode(p, r);
            return r;
        }

        vector<int64_t> decrypt(const Ciphertext& c) {
            Ciphertext c2 = c;
            if (c.isNttForm()) {
                evaluator->transformFromNttInplace(c2);
            }
            Plaintext p;
            decryptor->decrypt(c2, p);
            return decode(p);
        }

        void testSaveLoad(int repeatCount = 1000) {
            auto c1 = randomCiphertext();
            auto t1 = tim.registerTimer("Save Ciphertext");
            auto t2 = tim.registerTimer("Load Ciphertext");
            Ciphertext c2;
            for (int t = 0; t < repeatCount; t++) {
                std::stringstream ss;
                tim.tick(t1);
                c1.save(ss);
                tim.tock(t1);
                tim.tick(t2);
                c2.load(ss);
                tim.tock(t2);
            }
            auto p3 = randomPlaintext();
            auto t3 = tim.registerTimer("Save Plaintext");
            auto t4 = tim.registerTimer("Load Plaintext");
            Plaintext p4;
            for (int t = 0; t < repeatCount; t++) {
                std::stringstream ss;
                tim.tick(t3);
                p3.save(ss);
                tim.tock(t3);
                tim.tick(t4);
                p4.load(ss);
                tim.tock(t4);
            }
            printTimer(tim.gather(repeatCount));

            vector<int64_t> m1, m2, m3, m4;
            m1 = decrypt(c1);
            m2 = decrypt(c2);
            m3 = decode(p3);
            m4 = decode(p4);
            bool pass = is_equal(m1, m2);
            std::cout << pass_str(pass) << " | Save/Load Ciphertext" << std::endl;
            pass = is_equal(m3, m4);
            std::cout << pass_str(pass) << " | Save/Load Plaintext" << std::endl;
        }

        void testEncode(int repeatCount = 1000) {
            auto m1 = randomVector(slotCount, dataBound);
            // print m1
            std::cout << "m1: size=" << m1.size() << ", data="; printVector(m1);
            vector<int64_t> m2;
            auto t1 = tim.registerTimer("Encode");
            auto t2 = tim.registerTimer("Decode");
            for (int t = 0; t < repeatCount; t++) {
                Plaintext p1;
                tim.tick(t1);
                encoder->encode(m1, p1);
                tim.tock(t1);
                tim.tick(t2);
                encoder->decode(p1, m2);
                tim.tock(t2);
            }
            printTimer(tim.gather(repeatCount));

            // print m2
            std::cout << "m2: size=" << m2.size() << ", data="; printVector(m2);

            bool pass = is_equal(m1, m2);
            std::cout << pass_str(pass) << " | Encode/Decode" << std::endl;
        }
  
        void testAdd(int repeatCount = 1000) {
            auto c1 = randomCiphertext();
            auto c2 = randomCiphertext();
            Ciphertext c3;
            auto t1 = tim.registerTimer("Add-assign");
            auto t2 = tim.registerTimer("Add-inplace");
            for (int t = 0; t < repeatCount; t++) {
                tim.tick(t1);
                evaluator->add(c1, c2, c3);
                tim.tock(t1);
                tim.tick(t2);
                evaluator->addInplace(c3, c1);
                tim.tock(t2);
            }
            printTimer(tim.gather(repeatCount));

            vector<int64_t> m1, m2, m3;
            m1 = decrypt(c1);
            m2 = decrypt(c2);
            m3 = decrypt(c3);

            bool pass = is_equal(m3, vectorAdd(vectorAdd(m1, m2), m1));
            std::cout << pass_str(pass) << " | Add" << std::endl;
        }

        void testAddPlain(int repeatCount = 1000) {
            auto c1 = randomCiphertext();
            auto p2 = randomPlaintext();
            Ciphertext c3;
            auto t1 = tim.registerTimer("AddPlain-assign");
            auto t2 = tim.registerTimer("AddPlain-inplace");
            for (int t = 0; t < repeatCount; t++) {
                tim.tick(t1);
                evaluator->addPlain(c1, p2, c3);
                tim.tock(t1);
                tim.tick(t2);
                evaluator->addPlainInplace(c3, p2);
                tim.tock(t2);
            }
            printTimer(tim.gather(repeatCount));

            vector<int64_t> m1, m2, m3;
            m1 = decrypt(c1);
            m2 = decode(p2);
            m3 = decrypt(c3);

            bool pass = is_equal(m3, vectorAdd(vectorAdd(m1, m2), m2));
            std::cout << pass_str(pass) << " | AddPlain" << std::endl;
        }

        void testMultiplyPlain(int repeatCount = 1000) {
            auto c1 = randomCiphertext();
            auto p2 = randomPlaintext();
            Ciphertext c3;
            auto t1 = tim.registerTimer("MultiplyPlain-assign");
            auto t2 = tim.registerTimer("MultiplyPlain-inplace");
            for (int t = 0; t < repeatCount; t++) {
                tim.tick(t1);
                evaluator->multiplyPlain(c1, p2, c3);
                tim.tock(t1);
                tim.tick(t2);
                evaluator->multiplyPlainInplace(c3, p2);
                tim.tock(t2);
            }
            printTimer(tim.gather(repeatCount));

            vector<int64_t> m1, m2, m3;
            m1 = decrypt(c1);
            m2 = decode(p2);
            m3 = decrypt(c3);

            bool pass = is_equal(m3, vectorMultiply(vectorMultiply(m1, m2), m2));
            std::cout << pass_str(pass) << " | MultiplyPlain" << std::endl;
        }

        void testSquare(int repeatCount = 1000) {
            auto c1 = randomCiphertext();
            Ciphertext c2;
            Ciphertext c3;
            auto t1 = tim.registerTimer("Square-assign");
            auto t2 = tim.registerTimer("Square-inplace");
            for (int t = 0; t < repeatCount; t++) {
                tim.tick(t1);
                evaluator->square(c1, c2);
                tim.tock(t1);
                c3 = c1;
                tim.tick(t2);
                evaluator->squareInplace(c3);
                tim.tock(t2);
            }
            printTimer(tim.gather(repeatCount));

            vector<int64_t> m1, m2, m3;
            m1 = decrypt(c1);
            m2 = decrypt(c2);
            m3 = decrypt(c3);

            auto m1_2 = vectorMultiply(m1, m1);
            bool pass = is_equal(m2, m1_2);
            std::cout << pass_str(pass) << " | Square" << std::endl;
            pass = is_equal(m3, m1_2);
            std::cout << pass_str(pass) << " | SquareInplace" << std::endl;
        }


        void testMultiplyRescale(int repeatCount = 100) {
            auto c1 = randomCiphertext();
            auto c2 = randomCiphertext();
            Ciphertext c3, c4;
            Ciphertext c5;
            auto t1 = tim.registerTimer("Multiply-assign");
            auto t2 = tim.registerTimer("ModSwitch-assign");
            auto t2_ = tim.registerTimer("Relinearize-assign");
            auto t3 = tim.registerTimer("Multiply-inplace");
            auto t4 = tim.registerTimer("ModSwitch-inplace");
            auto t4_ = tim.registerTimer("Relinearize-inplace");
            for (int t = 0; t < repeatCount; t++) {
                tim.tick(t1);
                evaluator->multiply(c1, c2, c3);
                tim.tock(t1);
                tim.tick(t2);
                evaluator->modSwitchToNext(c3, c4);
                tim.tock(t2);
                c5 = c1;
                tim.tick(t3);
                evaluator->multiplyInplace(c5, c2);
                tim.tock(t3);
                tim.tick(t4);
                evaluator->modSwitchToNextInplace(c5);
                tim.tock(t4);
            }
            printTimer(tim.gather(repeatCount));

            vector<int64_t> m1, m2, m3, m4, m5;
            m1 = decrypt(c1);
            m2 = decrypt(c2);
            m3 = decrypt(c3);
            m4 = decrypt(c4);
            m5 = decrypt(c5);

            bool pass = is_equal(m3, vectorMultiply(m1, m2));
            std::cout << pass_str(pass) << " | Multiply" << std::endl;
            pass = is_equal(m4, m3);
            std::cout << pass_str(pass) << " | ModSwitch" << std::endl;
            pass = is_equal(m5, m3);
            std::cout << pass_str(pass) << " | Inplace" << std::endl;
        }

        void testRotateVector(int repeatCount = 100) {
            auto c1 = randomCiphertext();
            Ciphertext c2;
            auto t1 = tim.registerTimer("RotateRows-assign");
            auto t2 = tim.registerTimer("RotateRows-inplace");
            for (int t = 0; t < repeatCount; t++) {
                tim.tick(t1);
                evaluator->rotateRows(c1, 1, gk, c2);
                tim.tock(t1);
                tim.tick(t2);
                evaluator->rotateRowsInplace(c1, 1, gk);
                tim.tock(t2);
            }
            printTimer(tim.gather(repeatCount));

            vector<int64_t> m1, m2;
            m1 = decrypt(c1);
            m2 = decrypt(c2);
            bool pass = is_equal(m1, m2);
            std::cout << pass_str(pass) << " | rotateRows" << std::endl;
        }
        
        /*
            Test the performance of transformToNtt and transformToNttInplace
            args:
                repeatCount: the number of times to repeat the test
        */
        void testToNtt(int repeatCount = 1000) {
            auto c1 = randomCiphertext();
            Ciphertext c2;
            Ciphertext c3;
            auto t1 = tim.registerTimer("ToNtt-assign");
            auto t2 = tim.registerTimer("ToNtt-inplace");
            for (int t = 0; t < repeatCount; t++) {
                tim.tick(t1);
                // evaluator->square(c1, c2);
                evaluator->transformToNtt(c1, c2);
                tim.tock(t1);
                c3 = c1;
                tim.tick(t2);
                evaluator->transformToNttInplace(c3);
                tim.tock(t2);
            }
            printTimer(tim.gather(repeatCount));

            vector<int64_t> m1, m2, m3;
            m1 = decrypt(c1);
            m2 = decrypt(c2);
            m3 = decrypt(c3);
            bool pass = is_equal(m2, m1);
            std::cout << pass_str(pass) << " | ToNTT" << std::endl;
            pass = is_equal(m3, m1);
            std::cout << pass_str(pass) << " | ToNTTInplace" << std::endl;
        }

        /*
            Test the performance of transformFromNtt and transformFromNttInplace
            args:
                repeatCount: the number of times to repeat the test
        */
        void testFromNtt(int repeatCount = 1000) {
            auto c1 = randomCiphertext();
            evaluator->transformToNttInplace(c1);
            Ciphertext c2;
            Ciphertext c3;
            auto t1 = tim.registerTimer("FromNTT-assign");
            auto t2 = tim.registerTimer("FromNTT-inplace");
            for (int t = 0; t < repeatCount; t++) {
                tim.tick(t1);
                // evaluator->square(c1, c2);
                evaluator->transformFromNtt(c1, c2);
                tim.tock(t1);
                c3 = c1;
                tim.tick(t2);
                evaluator->transformFromNttInplace(c3);
                tim.tock(t2);
            }
            printTimer(tim.gather(repeatCount));

            vector<int64_t> m1, m2, m3;
            m1 = decrypt(c1);
            m2 = decrypt(c2);
            m3 = decrypt(c3);
            bool pass = is_equal(m2, m1);
            std::cout << pass_str(pass) << " | FromNTT" << std::endl;
            pass = is_equal(m3, m1);
            std::cout << pass_str(pass) << " | FromNTTInplace" << std::endl;
        }

        void testMultiplyPlainNtt(int repeatCount = 1000) {
            auto c1 = randomCiphertext();
            auto p2 = randomPlaintext();
            Plaintext p2_ntt;
            Ciphertext c3, c4, c5;
            auto t0 = tim.registerTimer("ToNttPlain-assign");
            auto t1 = tim.registerTimer("ToNtt-assign");
            auto t2 = tim.registerTimer("MultiplyPlainNTT-assign");
            auto t3 = tim.registerTimer("FromNtt-assign");
            auto t4 = tim.registerTimer("ToNttPlain-inplace");
            auto t5 = tim.registerTimer("ToNtt-inplace");
            auto t6 = tim.registerTimer("MultiplyPlainNTT-inplace");
            auto t7 = tim.registerTimer("FromNtt-inplace");
            for (int t = 0; t < repeatCount; t++) {
                tim.tick(t0);
                evaluator->transformToNtt(p2, c1.parmsID(), p2_ntt);
                tim.tock(t0);
                tim.tick(t1);
                evaluator->transformToNtt(c1, c3);
                tim.tock(t1);
                tim.tick(t2);
                evaluator->multiplyPlain(c3, p2_ntt, c4);
                tim.tock(t2);
                tim.tick(t3);
                evaluator->transformFromNtt(c4, c5);
                tim.tock(t3);
                c5 = c1;
                p2_ntt = p2;
                tim.tick(t4);
                evaluator->transformToNttInplace(p2_ntt, c1.parmsID());
                tim.tock(t4);
                tim.tick(t5);
                evaluator->transformToNttInplace(c5);
                tim.tock(t5);
                tim.tick(t6);
                evaluator->multiplyPlainInplace(c5, p2_ntt);
                tim.tock(t6);
                tim.tick(t7);
                evaluator->transformFromNttInplace(c5);
                tim.tock(t7);
            }
            printTimer(tim.gather(repeatCount));

            vector<int64_t> m1, m2, m5;
            m1 = decrypt(c1);
            m2 = decode(p2);
            m5 = decrypt(c5);
            bool pass = is_equal(m5, vectorMultiply(m1, m2));
            std::cout << pass_str(pass) << " | MultiplyPlainNTT" << std::endl;
        }

        void testAll() {
            this->testSaveLoad();
            this->testEncode();
            this->testEncrypt();
            this->testAdd();
            this->testAddPlain();
            this->testMultiplyRescale();
            this->testMultiplyPlain();
            this->testSquare();
            this->testRotateVector();
            this->testToNtt();
            this->testFromNtt();
            this->testMultiplyPlainNtt();
            this->testMemoryPool();
        }
    };
}

int main() {

    std::vector<int> qs = {60, 60, 60};
    int plainModBits = 20;
    int polyModDeg = 8192;

    std::cout << "----- TimeTest cuda CKKS -----\n";
    troytest::TimeTestCKKS test(polyModDeg, qs);
    test.testAll();

    std::cout << "----- TimeTest cuda BFV -----\n";
    troytest::TimeTestBFVBGV test2(false, polyModDeg, plainModBits, qs);
    test2.testAll();

    std::cout << "----- TimeTest cuda BGV -----\n";
    troytest::TimeTestBFVBGV test3(true, polyModDeg, plainModBits, qs);
    test3.testAll();
    return 0;
}