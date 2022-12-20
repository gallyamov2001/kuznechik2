#include <chrono>
#include <iostream>
#include <string>
#include "reader.hpp"


int main()
{
    Mode mode = Mode::ECB;
    std::string file("input.txt");

    Encryptor E(mode);

    // encryption
    uint64_t num_blocks = E.ReadText(file);

    Timer timer;
    // encrypt
    timer.Start();

    E.Encrypt();

    timer.Finish();
    std::cout << "Encrypt time: " << timer.GetMilliseconds() << " ms\n";
    std::cout << "Speed: " << num_blocks * 1000.0  / (timer.GetMilliseconds() * 65536) << "MB/s\n";

    E.SaveText("encrypted.txt");

    // decryption
    E.ReadText("encrypted.txt", true);

    timer.Start();

    E.Decrypt();

    timer.Finish();
    std::cout << "Decrypt time: " << timer.GetMilliseconds() << " ms\n";
    std::cout << "Speed: " << num_blocks * 1000.0  / (timer.GetMilliseconds() * 65536) << "MB/s\n";

    E.SaveText("decrypted.txt", true);

    return 0;
}
