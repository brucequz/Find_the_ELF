# Find_the_ELF


Hi! Currently this readme file does not serve the purpose to be a descriptive illustration of this on-going project. It's just mainly a space where I put down my thoughts and meeting notes.

I do not imagine anyone else would find this unorganized text useful. But have a nice day!


1. crc convolution and deconvolution have replaced polynomial multiplication. They seem to be the same on paper but different in code. Might investigate later.
2. still using checkCRC (if all remainder are zero) for both ZT and TB list decoder. Worked so far.
3. As of Oct.4th, the main source code is inside test/DuoCodec_Test. main.cpp is outdated.
4. As of Oct.4th, I am considering migrating source code to main.cpp and getting ready for Monte Carlo simulation.
5. As of Oct.21th, I have refactored the Dual List Decoder and implemented dual List Map data structure to achieve O(1) search. std::unordered_map cannot hash std::vector<int> as keys so I had to use std::map.
6. As of Oct.21th, softViterbiDecoding for 2^14 states is too slow. Running 1e4 experiments might take overnight. Had to verify this with beryl.



N states, T stage
 Time : O(N^2 * T)


ideas for the structure
I am going to work with Message Information directly, or with a vector of messageInformation.



