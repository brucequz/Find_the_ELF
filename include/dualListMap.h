#ifndef DUALLISTMAP_H
#define DUALLISTMAP_H

#include <map>
#include <queue>
#include <vector>
#include <iostream>

struct MessageInformation;

struct DLDInfo {
  double combined_metric;
  std::vector<int> message;
  std::vector<int> list_ranks;
  std::vector<double> received_signal;
};

// Define a custom comparison function for the priority queue
struct CompareCombinedMetric {
  bool operator()(const DLDInfo& a, const DLDInfo& b) const {
    return a.combined_metric >
           b.combined_metric;  // Lower combined_metric at the top
  }
};

class DualListMap{
  public:
    DualListMap() {};
    ~DualListMap() {};

    void insert(const MessageInformation& mi);  // insert to priority queue when agreed message is found
    int queue_size() {return agreed_messages_.size();};
    DLDInfo pop_queue();


  private:
    std::map<std::vector<int>, MessageInformation> dual_list_map_;
    std::priority_queue<DLDInfo, std::vector<DLDInfo>, CompareCombinedMetric> agreed_messages_;

};


#endif


