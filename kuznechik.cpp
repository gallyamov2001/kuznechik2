#include <algorithm>
#include <functional>
#include <cstdio>
#include <cassert>
#include <xmmintrin.h>
#include <immintrin.h>
#include <emmintrin.h>
#include "kuznechik.hpp"
#include "tables.hpp"

using namespace std;

Kuznechik::Kuznechik()
{
    MakeMulTable();
    MakeCoefTable();
    MakeEncTable();
    MakeDecTable();
}

__m128i& M128B(uint8_t *ptr)
{
    return *reinterpret_cast<__m128i*>(ptr);
}

void Kuznechik::Encrypt(Data& data, const Key& key)
{
    KeyPair current_key;
    copy(key.begin(), key.begin() + block_size, current_key.first.begin());
    copy(key.begin() + block_size, key.end(), current_key.second.begin());

    Keys keys = MakeKeys(current_key);
    Block feedback = current_key.first;

    for (auto& block: data)
    {
        EncryptBlock(feedback, keys);
        FuncX(block, feedback);
        feedback = block;
    }
}

void Kuznechik::Decrypt(Data& data, const Key& key)
{
    KeyPair current_key;
    copy(key.begin(), key.begin() + block_size, current_key.first.begin());
    copy(key.begin() + block_size, key.end(), current_key.second.begin());

    Keys keys = MakeKeys(current_key);

    Block feedback = current_key.first;

    for (size_t i = 0; i + 1 < data.size(); i += 2)
    {
        DoubleBlock d_block = _mm256_castsi128_si256(feedback);
        d_block = _mm256_inserti128_si256(d_block, data[i], 1);
        EncryptBlock(d_block, keys);
        feedback = data[i + 1];
        __m128i tmp1 = _mm256_extracti128_si256 (d_block, 0);
        __m128i tmp2 = _mm256_extracti128_si256 (d_block, 1);
        FuncX(data[i]    , tmp1);
        FuncX(data[i + 1], tmp2);
    }
    if (data.size() % 2)
    {
        EncryptBlock(feedback, keys);
        FuncX(data.back(), feedback);
    }
    
}

static Kuznechik::Block mask =
{{
    0x00, 0xFF, 0x00, 0xFF, 0x00, 0xFF, 0x00, 0xFF,
    0x00, 0xFF, 0x00, 0xFF, 0x00, 0xFF, 0x00, 0xFF
}};

static Kuznechik::DoubleBlock d_mask =
{{
    0x00, 0xFF, 0x00, 0xFF, 0x00, 0xFF, 0x00, 0xFF,
    0x00, 0xFF, 0x00, 0xFF, 0x00, 0xFF, 0x00, 0xFF,
    0x00, 0xFF, 0x00, 0xFF, 0x00, 0xFF, 0x00, 0xFF,
    0x00, 0xFF, 0x00, 0xFF, 0x00, 0xFF, 0x00, 0xFF

}};


template <typename BlockType>
void Kuznechik::EncryptBlock(BlockType& data, const Keys& keys)
{

#ifdef VERBOSE
    printf("keys:\n");
    for (const auto& key: keys)
        key.Dump();
    printf("before encryption:\n");
    data.Dump();
#endif

    for (unsigned i = 0; i < num_rounds - 1; ++i)
        FuncXSL(data, keys[i]);
    FuncX(data, keys[num_rounds - 1]);

#ifdef VERBOSE
    printf("after encryption:\n");
    data.Dump();
#endif
}

template <typename BlockType>
void Kuznechik::DecryptBlock(BlockType& data, const Keys& keys)
{

#ifdef VERBOSE
    printf("keys:\n");
    for (const auto& key: keys)
        key.Dump();
    printf("before decryption:\n");
    data.Dump();
#endif

    transform(data.begin(), data.end(), data.begin(),
                   [](uint8_t idx) { return Table::S[idx]; });
    for (unsigned i = num_rounds - 1; i > 0; --i)
        FuncInvXLS(data, keys[i]);
    transform(data.begin(), data.end(), data.begin(),
                   [](uint8_t idx) { return Table::invS[idx]; });
    FuncX(data, keys[0]);

#ifdef VERBOSE
    printf("after decryption:\n");
    data.Dump();
#endif
}

template<typename BlockType>
void Kuznechik::FuncXSL(BlockType& data, const Block& key)
{
    FuncX(data, key);
    FuncLS(data, enc_ls_table);
}

template<typename BlockType>
void Kuznechik::FuncInvXLS(BlockType& data, const Block& key)
{
    FuncLS(data, dec_ls_table);
    FuncX(data, key);
}

__attribute__((always_inline))
void Kuznechik::FuncLS(Block& data, const LookupTable& lut)
{
    __m128i tmp1 = _mm_and_si128 (mask, data);
    __m128i tmp2 = _mm_andnot_si128 (mask, data);
    tmp1 = _mm_srli_epi64 (tmp1, 4);
    tmp2 = _mm_slli_epi64 (tmp2, 4);

    uint8_t* table = (uint8_t*)&lut;
    __m128i vector1 = _mm_load_si128(reinterpret_cast<const __m128i*>(table + _mm_extract_epi16(tmp2, 0) + 0x0000));
    __m128i vector2 = _mm_load_si128(reinterpret_cast<const __m128i*>(table + _mm_extract_epi16(tmp1, 0) + 0x1000));

    vector1 = _mm_xor_si128(vector1, M128B(table + _mm_extract_epi16(tmp2, 1) + 0x2000));
    vector2 = _mm_xor_si128(vector2, M128B(table + _mm_extract_epi16(tmp1, 1) + 0x3000));

    vector1 = _mm_xor_si128(vector1, M128B(table + _mm_extract_epi16(tmp2, 2) + 0x4000));
    vector2 = _mm_xor_si128(vector2, M128B(table + _mm_extract_epi16(tmp1, 2) + 0x5000));

    vector1 = _mm_xor_si128(vector1, M128B(table + _mm_extract_epi16(tmp2, 3) + 0x6000));
    vector2 = _mm_xor_si128(vector2, M128B(table + _mm_extract_epi16(tmp1, 3) + 0x7000));

    vector1 = _mm_xor_si128(vector1, M128B(table + _mm_extract_epi16(tmp2, 4) + 0x8000));
    vector2 = _mm_xor_si128(vector2, M128B(table + _mm_extract_epi16(tmp1, 4) + 0x9000));

    vector1 = _mm_xor_si128(vector1, M128B(table + _mm_extract_epi16(tmp2, 5) + 0xA000));
    vector2 = _mm_xor_si128(vector2, M128B(table + _mm_extract_epi16(tmp1, 5) + 0xB000));

    vector1 = _mm_xor_si128(vector1, M128B(table + _mm_extract_epi16(tmp2, 6) + 0xC000));
    vector2 = _mm_xor_si128(vector2, M128B(table + _mm_extract_epi16(tmp1, 6) + 0xD000));

    vector1 = _mm_xor_si128(vector1, M128B(table + _mm_extract_epi16(tmp2, 7) + 0xE000));
    vector2 = _mm_xor_si128(vector2, M128B(table + _mm_extract_epi16(tmp1, 7) + 0xF000));
    data = _mm_xor_si128(vector1, vector2);
}

__attribute__((always_inline))
void Kuznechik::FuncLS(DoubleBlock& data, const LookupTable& lut)
{
    __m256i tmp1 = _mm256_and_si256 (d_mask, data);
    __m256i tmp2 = _mm256_andnot_si256 (d_mask, data);

    tmp1 = _mm256_srli_epi64(tmp1, 4);
    tmp2 = _mm256_slli_epi64(tmp2, 4);

    __m256i vector1 = _mm256_setzero_si256 ();
    __m256i vector2 = _mm256_setzero_si256 ();
    __m256i vector3 = _mm256_setzero_si256 ();
    data = _mm256_setzero_si256 ();
    uint8_t *table = (uint8_t*)&lut;

    data = _mm256_inserti128_si256(data, M128B(table + _mm256_extract_epi16(tmp2, 0) + 0x0000), 0);
    data = _mm256_inserti128_si256(data, M128B(table + _mm256_extract_epi16(tmp2, 8) + 0x0000), 1);

    vector2 = _mm256_inserti128_si256(vector2, M128B(table + _mm256_extract_epi16(tmp2, 1) + 0x2000), 0);
    vector2 = _mm256_inserti128_si256(vector2, M128B(table + _mm256_extract_epi16(tmp2, 9) + 0x2000), 1);

    vector1 = _mm256_inserti128_si256(vector1, M128B(table + _mm256_extract_epi16(tmp1, 0) + 0x1000), 0);
    vector1 = _mm256_inserti128_si256(vector1, M128B(table + _mm256_extract_epi16(tmp1, 8) + 0x1000), 1);
    
    vector3 = _mm256_inserti128_si256(vector3, M128B(table + _mm256_extract_epi16(tmp1, 1) + 0x3000), 0);
    vector3 = _mm256_inserti128_si256(vector3, M128B(table + _mm256_extract_epi16(tmp1, 9) + 0x3000), 1);

    data = _mm256_xor_si256(data, vector2);
    vector1 = _mm256_xor_si256(vector1, vector3);

    vector2 = _mm256_inserti128_si256(vector2, M128B(table + _mm256_extract_epi16(tmp2, 2) + 0x4000), 0);
    vector2 = _mm256_inserti128_si256(vector2, M128B(table + _mm256_extract_epi16(tmp2,10) + 0x4000), 1);

    vector3 = _mm256_inserti128_si256(vector3, M128B(table + _mm256_extract_epi16(tmp1, 2) + 0x5000), 0);
    vector3 = _mm256_inserti128_si256(vector3, M128B(table + _mm256_extract_epi16(tmp1, 10) + 0x5000), 1);

    data = _mm256_xor_si256(data, vector2);
    vector1 = _mm256_xor_si256(vector1, vector3);

    vector2 = _mm256_inserti128_si256(vector2, M128B(table + _mm256_extract_epi16(tmp2, 3) + 0x6000), 0);
    vector2 = _mm256_inserti128_si256(vector2, M128B(table + _mm256_extract_epi16(tmp2,11) + 0x6000), 1);

    vector3 = _mm256_inserti128_si256(vector3, M128B(table + _mm256_extract_epi16(tmp1, 3) + 0x7000), 0);
    vector3 = _mm256_inserti128_si256(vector3, M128B(table + _mm256_extract_epi16(tmp1, 11) + 0x7000), 1);

    data = _mm256_xor_si256(data, vector2);
    vector1 = _mm256_xor_si256(vector1, vector3);

    vector2 = _mm256_inserti128_si256(vector2, M128B(table + _mm256_extract_epi16(tmp2, 4) + 0x8000), 0);
    vector2 = _mm256_inserti128_si256(vector2, M128B(table + _mm256_extract_epi16(tmp2,12) + 0x8000), 1);

    vector3 = _mm256_inserti128_si256(vector3, M128B(table + _mm256_extract_epi16(tmp1, 4) + 0x9000), 0);
    vector3 = _mm256_inserti128_si256(vector3, M128B(table + _mm256_extract_epi16(tmp1, 12) + 0x9000), 1);

    data = _mm256_xor_si256(data, vector2);
    vector1 = _mm256_xor_si256(vector1, vector3);

    vector2 = _mm256_inserti128_si256(vector2, M128B(table + _mm256_extract_epi16(tmp2, 5) + 0xA000), 0);
    vector2 = _mm256_inserti128_si256(vector2, M128B(table + _mm256_extract_epi16(tmp2,13) + 0xA000), 1);

    vector3 = _mm256_inserti128_si256(vector3, M128B(table + _mm256_extract_epi16(tmp1, 5) + 0xB000), 0);
    vector3 = _mm256_inserti128_si256(vector3, M128B(table + _mm256_extract_epi16(tmp1, 13) + 0xB000), 1);

    data = _mm256_xor_si256(data, vector2);
    vector1 = _mm256_xor_si256(vector1, vector3);

    vector2 = _mm256_inserti128_si256(vector2, M128B(table + _mm256_extract_epi16(tmp2, 6) + 0xC000), 0);
    vector2 = _mm256_inserti128_si256(vector2, M128B(table + _mm256_extract_epi16(tmp2,14) + 0xC000), 1);

    vector3 = _mm256_inserti128_si256(vector3, M128B(table + _mm256_extract_epi16(tmp1, 6) + 0xD000), 0);
    vector3 = _mm256_inserti128_si256(vector3, M128B(table + _mm256_extract_epi16(tmp1, 14) + 0xD000), 1);

    data = _mm256_xor_si256(data, vector2);
    vector1 = _mm256_xor_si256(vector1, vector3);

    vector2 = _mm256_inserti128_si256(vector2, M128B(table + _mm256_extract_epi16(tmp2, 7) + 0xE000), 0);
    vector2 = _mm256_inserti128_si256(vector2, M128B(table + _mm256_extract_epi16(tmp2,15) + 0xE000), 1);

    vector3 = _mm256_inserti128_si256(vector3, M128B(table + _mm256_extract_epi16(tmp1, 7) + 0xF000), 0);
    vector3 = _mm256_inserti128_si256(vector3, M128B(table + _mm256_extract_epi16(tmp1, 15) + 0xF000), 1);

    data = _mm256_xor_si256(data, vector2);
    vector1 = _mm256_xor_si256(vector1, vector3);
    data = _mm256_xor_si256(data, vector1);
}

Kuznechik::Keys Kuznechik::MakeKeys(const KeyPair& key)
{
    Keys keys{{key.first, key.second}};
    for (unsigned i = 1; i < num_rounds / 2; ++i)
    {
        keys[2 * i] = keys[2 * i - 2];
        keys[2 * i + 1] = keys[2 * i - 1];
        for (unsigned j = 0; j < 8; ++j)
            FuncF(keys[2 * i], keys[2 * i + 1], coef_table[i - 1][j]);
    }
    return keys;
}

void Kuznechik::FuncF(Block& data1, Block& data0, const Block& key)
{
    Block tmp = data1;
    FuncXSL(tmp, key);
    FuncX(tmp, data0);
    tie(data1, data0) = make_pair(tmp, data1);
}

void Kuznechik::MakeMulTable()
{
    for (unsigned i = 0; i < 256; ++i)
        for (unsigned j = 0; j < 256; ++j)
            mul_table[i][j] = PolyMul(i, j);
}

void Kuznechik::MakeCoefTable()
{
    for (unsigned i = 0; i < num_rounds / 2 - 1; ++i)
        for (unsigned j = 0; j < 8; ++j)
        {
            Block& C = coef_table[i][j];
            C.fill(0);
            C.back() = i * 8 + j + 1;
            FuncL(C);
        }
}

void Kuznechik::MakeEncTable()
{
    //initializing the martix
    Matrix l_matrix;
    for (size_t i = 0; i < block_size; ++i)
        for (size_t j = 0; j < block_size; ++j)
            if (i == 0)
                l_matrix[i][j] = Table::lin[j];
            else if (i == j + 1)
                l_matrix[i][j] = 1;
            else
                l_matrix[i][j] = 0;
    // power l_matrix to degree of 2^4
    for (unsigned i = 0; i < 4; ++i)
        l_matrix = SqrMatrix(l_matrix);

    for (size_t i = 0; i < block_size; ++i)
        for (size_t j = 0; j < 256; ++j)
            for (size_t k = 0; k < block_size; ++k)
                enc_ls_table[i][j][k] = mul_table[Table::S[j]][l_matrix[k][i]];
}

void Kuznechik::MakeDecTable()
{
    //initializing the martix
    Matrix l_matrix;
    for (size_t i = 0; i < block_size; ++i)
        for (size_t j = 0; j < block_size; ++j)
            if (i == block_size - 1)
                l_matrix[i][j] = Table::lin[(j + 15) % block_size];
            else if (i + 1 == j)
                l_matrix[i][j] = 1;
            else
                l_matrix[i][j] = 0;
    // power l_matrix to degree of 2^4
    for (unsigned i = 0; i < 4; ++i)
        l_matrix = SqrMatrix(l_matrix);

    for (size_t i = 0; i < block_size; ++i)
        for (size_t j = 0; j < 256; ++j)
            for (size_t k = 0; k < block_size; ++k)
                dec_ls_table[i][j][k] = mul_table[Table::invS[j]][l_matrix[k][i]];
}

uint8_t Kuznechik::PolyMul (uint8_t left, uint8_t right)
{
    // p(x) = x^8 + x^7 + x^6 + x + 1 => 0b111000011 => 0xC3 (without MSB)
    uint8_t res = 0;
    while (left && right)
    {
        if (right & 1)
            res ^= left;
        left = (left << 1) ^ (left & 0x80 ? 0xC3 : 0x00);
        right >>= 1;
    }
    return res;
}

Kuznechik::Matrix Kuznechik::SqrMatrix(const Matrix& mat)
{
    Matrix res{};
    for (size_t i = 0; i < block_size; ++i)
        for (size_t j = 0; j < block_size; ++j)
            for (size_t k = 0; k < block_size; ++k)
                res[i][j] ^= mul_table[mat[i][k]][mat[k][j]];
    return res;
}

void Kuznechik::FuncX(Block& data, const Block& key)
{
    data = _mm_xor_si128(data, key);
}

void Kuznechik::FuncX(DoubleBlock& data, const Block& key)
{
    __m256i d_key = _mm256_castsi128_si256(key);
    d_key = _mm256_inserti128_si256 (d_key, key, 1);
    data = _mm256_xor_si256(data, d_key);
}

void Kuznechik::FuncL(Block& data)
{
    for (unsigned i = 0; i < block_size; ++i)
    {
        uint8_t tmp = 0;
        for (size_t j = 0; j < block_size; ++j)
            tmp ^= mul_table[data[j]][Table::lin[j]];
        copy_backward(data.begin(), data.end() - 1, data.end());
        data[0] = tmp;
    }
}

void Kuznechik::FuncInvL(Block& data)
{
    for (unsigned i = 0; i < block_size; ++i)
    {
        uint8_t elem = data[0];
        copy(data.begin() + 1, data.end(), data.begin());
        data[block_size - 1] = elem;

        uint8_t tmp = 0;
        for (size_t j = 0; j < block_size; ++j)
            tmp ^= mul_table[data[j]][Table::lin[j]];
        data[15] = tmp;
    }
}
