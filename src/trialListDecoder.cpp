#include "../include/trialListDecoder.h"

#include "../include/viterbiDecoder.h"
#include "../include/feedForwardTrellis.h"
#include "../include/CONSTANTS.h"

// void TrialListMap::insert(const MessageInformation& mi) {
//   auto it = trial_list_map_.find(mi.message);
//   if (it == trial_list_map_.end()) {
//     if (mi.decoder_index == -1) {
//       throw std::invalid_argument("Decoder index not set");
//     }
//     TLDInfo new_tld_info;
//     new_tld_info.mi = mi;
//     new_tld_info.decoder_index[mi.decoder_index] = true;
//     trial_list_map_.insert({mi.message, new_tld_info});
//   } else {
//     // key already exists, update the decoder index
//     it->second.decoder_index[mi.decoder_index] = true;
//     // three way match is found
//     if (it->second.decoder_index[0] && it->second.decoder_index[1] && it->second.decoder_index[2]) {
//       agreed_messages_.push(it->second);
//       trial_list_map_.erase(it);
//     }
//   }
// }

// TLDInfo TrialListMap::pop_queue() {
//   TLDInfo top = agreed_messages_.top();
//   agreed_messages_.pop();
//   return top;
// }

// TLDInfo TrialListMap::get_top() {
//   return agreed_messages_.top();
// }

TrialListDecoder::TrialListDecoder(CodeInformation encoder, std::vector<CodeInformation> code_info, int max_path_to_search) 
: ViterbiDecoder(encoder, code_info, max_path_to_search) {
  if (code_info.size() != 3) {
    throw std::invalid_argument("TrialListDecoder only supports 3 codes");
  }
  for (auto& code : code_info_) {
    FeedForwardTrellis* trellis = new FeedForwardTrellis(code);
    trellis_ptrs_.push_back(trellis);
  }

}

TrialListDecoder::~TrialListDecoder() {
  for (auto& trellis : trellis_ptrs_) {
    delete trellis;
  }
  delete encoder_trellis_ptr_;
}

TLDInfo TrialListDecoder::TLD_BAM(std::vector<double> received_signal,  std::map<std::string, std::chrono::milliseconds>& timeDurations) {

  double best_current_match = INT_MAX;
  std::vector<double> received_codec_0;
  std::vector<double> received_codec_1;
  std::vector<double> received_codec_2;

  for (size_t i = 0; i < received_signal.size(); i += 3) {
    if (i % 3 == 0) {
      received_codec_0.push_back(received_signal[i]);
      received_codec_0.push_back(received_signal[i + 1]);
    } else if (i % 3 == 1) {
      received_codec_1.push_back(received_signal[i]);
      received_codec_1.push_back(received_signal[i + 1]);
    } else {
      received_codec_2.push_back(received_signal[i - 2]);
      received_codec_2.push_back(received_signal[i]); 
    }
  }
  assert(received_codec_0.size() == BLOCK_SIZE);
  assert(received_codec_1.size() == BLOCK_SIZE);
  assert(received_codec_2.size() == BLOCK_SIZE);

  bool best_combined_found = false;
  CodeInformation code_0 = code_info_[0];
  CodeInformation code_1 = code_info_[1];
  CodeInformation code_2 = code_info_[2];

  TLDInfo unmatched_best;

  return unmatched_best;
}
