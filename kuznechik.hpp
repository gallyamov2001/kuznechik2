#ifndef Kuznechik_H
#define Kuznechik_H

#include <array>
#include <vector>
#include <utility>
#include <cstdlib>
#include <cstdint>
#include <immintrin.h>

using namespace std;

template<typename T>
class AlignedAllocator
{
public:
    typedef T value_type;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;

    T* allocate(size_t n)
    {
        if (!n)
            return nullptr;
        if (n > max_size)
            throw length_error("AlignedAllocator::allocate() - integer overflow");
        void *p = aligned_alloc(sizeof(T) * 2, n * sizeof(T));
        if (!p)
            throw bad_alloc();
        return static_cast<T*>(p);
    }

    void deallocate(T* p, size_t n __attribute__((unused)))
    {
        free(static_cast<void*>(p));
    }

    bool operator ==(const AlignedAllocator& rhs __attribute__((unused))) const
    {
        return true;
    }

    bool operator !=(const AlignedAllocator& rhs) const
    {
        return !(*this == rhs);
    }

private:
    static constexpr size_t max_size =
        (static_cast<size_t>(0) - static_cast<size_t>(1)) / sizeof(T);
};

enum class Mode
{
    ECB
};

class Kuznechik
{
public:
    static constexpr size_t block_size = 16;
    static constexpr unsigned num_rounds = 10;

    struct alignas(block_size) Block : array<uint8_t, block_size>
    {
        Block() = default;

        Block(const __m128i& val)
        {
            *reinterpret_cast<__m128i*>(this) = val;
        }

        Block(const array<uint8_t, block_size>& arr)
            : array<uint8_t, block_size>(arr)
        {
        }

        operator __m128i() const
        {
            return *reinterpret_cast<const __m128i*>(this);
        }

        void Dump()
        {
            for (auto x: *this)
                printf("%02x", x);
            printf("\n");
        }
    };

    struct alignas(block_size * 2) DoubleBlock : array<uint8_t, block_size * 2>
    {
        DoubleBlock() = default;

        DoubleBlock(const __m256i& val)
        {
            *reinterpret_cast<__m256i*>(this) = val;
        }

        DoubleBlock(const array<uint8_t, block_size * 2>& arr)
            : array<uint8_t, block_size * 2>(arr)
        {
        }

        operator __m256i() const
        {
            return *reinterpret_cast<const __m256i*>(this);
        }

        void Dump()
        {
            for (auto x: *this)
                printf("%02x", x);
            printf("\n");
        }
    };

    using Key = array<uint8_t, block_size * 2>;
    using Data = vector<Block, AlignedAllocator<Block>>;

    Kuznechik();

    void Encrypt(Data& data, const Key& key);
    void Decrypt(Data& data, const Key& key);
private:
    using Matrix = array<Block, block_size>;
    using Keys = array<Block, num_rounds>;
    using KeyPair = pair<Block, Block>;
    using LookupTable = Block[block_size][256];
    Block coef_table[num_rounds / 2 - 1][8];
    uint8_t mul_table[256][256];
    LookupTable enc_ls_table;
    LookupTable dec_ls_table;
    
    void MakeMulTable();
    void MakeCoefTable();
    void MakeEncTable();
    void MakeDecTable();

    template<typename BlockType>
    void EncryptBlock(BlockType& data, const Keys& keys);
    template<typename BlockType>
    void DecryptBlock(BlockType& data, const Keys& keys);

    template<typename BlockType>
    void FuncXSL(BlockType& data, const Block& key);
    template<typename BlockType>
    void FuncInvXLS(BlockType& data, const Block& key);
    void FuncLS(Block& data, const LookupTable& lut);
    void FuncLS(DoubleBlock& data, const LookupTable& lut);

    Keys MakeKeys(const KeyPair& key);
    void FuncF(Block& data0, Block& data1, const Block& key);
    void FuncX(Block& data, const Block& key);
    void FuncX(DoubleBlock& data, const Block& key);
    void FuncL(Block& data);
    void FuncInvL(Block& data);


    uint8_t PolyMul(uint8_t left, uint8_t right);
    Matrix SqrMatrix(const Matrix& mat);
};

#endif
