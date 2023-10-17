#include <gtest/gtest.h>
#include "../include/minHeap.h"

TEST(MinHeapTest, HeapOperation) {
  MinHeap heap;
  DetourNode node1;
  node1.path_metric = 13.19048019070684;
  DetourNode node2;
  node2.path_metric = 12.03847712458178;
  DetourNode node3;
  node3.path_metric = 21.163543088833173;
  DetourNode node4;
  node4.path_metric = 15.290142686594368;
  DetourNode node5;
  node5.path_metric = 19.159659801394159;
  DetourNode node6;
  node6.path_metric = 1.3149828534360761;
  DetourNode node7;
  node7.path_metric = 17.556868205710984;
  DetourNode node8;
  node8.path_metric = 13.9989;

  heap.insert(node1);
  heap.insert(node2);
  heap.insert(node3);
  heap.insert(node4);
  heap.insert(node5);
  heap.insert(node6);
  heap.insert(node7);
  heap.insert(node8);

  EXPECT_EQ(heap.size(), 8);
  EXPECT_EQ(heap.pop().path_metric, 1.3149828534360761);
  //EXPECT_EQ(heap.size(), 7);
  EXPECT_EQ(heap.pop().path_metric, 12.03847712458178);
  //EXPECT_EQ(heap.size(), 6);
  EXPECT_EQ(heap.pop().path_metric, 13.19048019070684);
  //EXPECT_EQ(heap.size(), 5);
  EXPECT_EQ(heap.pop().path_metric, 13.9989);
  //EXPECT_EQ(heap.size(), 4);
  EXPECT_EQ(heap.pop().path_metric, 15.290142686594368);
  //EXPECT_EQ(heap.size(), 3);
  EXPECT_EQ(heap.pop().path_metric, 17.556868205710984); 
  //EXPECT_EQ(heap.size(), 2);
  EXPECT_EQ(heap.pop().path_metric, 19.159659801394159);

  EXPECT_EQ(heap.pop().path_metric, 21.163543088833173);

}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}