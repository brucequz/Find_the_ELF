#ifndef TRIALLISTDECODER_H
#define TRIALLISTDECODER_H

#include <iostream>
#include <vector>
#include <map>
#include <queue>

#include "viterbiDecoder.h"

struct FeedForwardTrellis;

namespace trial_list_decoder_utils {

}

struct TLDInfo {
  MessageInformation mi;
  std::vector<bool> decoder_index;
};



class TrialListMap {
  public:
    TrialListMap() {};
    ~TrialListMap() {};

    void insert(const MessageInformation& mi);
    int queue_size() {return agreed_messages_.size();};
    TLDInfo pop_queue();
    TLDInfo get_top();

  private:  
    std::map<std::vector<int>, TLDInfo> trial_list_map_; // dictionary
    std::priority_queue<TLDInfo, std::vector<TLDInfo>, CompareCombinedMetric<TLDInfo>> agreed_messages_; // priority queue
};


class TrialListDecoder : private ViterbiDecoder {
  public:
    TrialListDecoder(CodeInformation encoder, std::vector<CodeInformation> code_info, int max_path_to_search);
    ~TrialListDecoder();

    TLDInfo TLD_BAM(std::vector<double> received_signal, std::map<std::string, std::chrono::milliseconds>& timeDurations);

};

#endif