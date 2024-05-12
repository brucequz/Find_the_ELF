#include "../include/dualListMap.h"

#include "../include/dualListDecoder.h"

void DualListMap::insert(const MessageInformation& mi) {
  /*
  Insert a MessageInformation into an unordered map
  Insertion logic:

  extract the message from input argument
  search if there exists a key-value pair with
  if there is a match
    add the maximum likelihood MI to the priority queue
    erase that message from the unordered map
    if dual list decoder gets the maximum likelihood MI,
    set BCM and thresholds for list decoder 0 and 1.
  end if

  */
  auto it = dual_list_map_.find(mi.message); // finding the match in the dictionary
  if (it != dual_list_map_.end()) {
    // key exists
    // if there is a match
    DLDInfo agreed_message;
    agreed_message.combined_metric = mi.path_metric;// + it->second.path_metric;
    if (it->second.decoder_index == 0) {
      agreed_message.list_ranks = {it->second.list_rank, mi.list_rank};
      // std::cout << "DLD 0 result is already in there with " <<
      // it->second.path_metric << " pathMetric" << std::endl;
    } else if (it->second.decoder_index == 1) {
      agreed_message.list_ranks = {mi.list_rank, it->second.list_rank};
      // std::cout << "DLD 1 result is already in there with " <<
      // it->second.path_metric << " pathMetric" << std::endl;
    } else {
      std::cerr << "Invalid decoder order" << std::endl;
    }
    // std::cout << "The newly inserted one has " << mi.path_metric << "as
    // pathMetric." << std::endl;
    agreed_message.message = mi.message;
    agreed_message.symbol_metrics = mi.symbol_metrics;
    agreed_messages_.push(agreed_message);
    dual_list_map_.erase(mi.message);
  } else {
    // key does not exist
    dual_list_map_[mi.message] = mi;
  }
}

DLDInfo DualListMap::pop_queue() {
  if (queue_size() == 0) {
    std::cerr << "Invalid access to empty queue" << std::endl;
  }
  DLDInfo output = agreed_messages_.top();
  agreed_messages_.pop();
  return output;
}


DLDInfo DualListMap::get_top() {
  if (queue_size() == 0) {
    std::cerr << "Invalid access to empty queue" << std::endl;
  }
  DLDInfo output = agreed_messages_.top();
  return output;
}
