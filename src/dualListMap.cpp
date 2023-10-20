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
  auto it = dual_list_map_.find(mi.message);
  if (it != dual_list_map_.end()) {
    // key exists
    DLDInfo agreed_message;
    agreed_message.combined_metric = mi.path_metric + it->second.path_metric;
    // be cautious since the existing message could come from any list decoder
    if (it->second.decoder_index == 0) {
      agreed_message.list_ranks = {it->second.list_rank, mi.list_rank};
    } else if (it->second.decoder_index == 1) {
      agreed_message.list_ranks = {mi.list_rank, it->second.list_rank};
    } else {
      std::cerr << "Invalid decoder order" << std::endl;
    }
    agreed_message.message = mi.message;
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
