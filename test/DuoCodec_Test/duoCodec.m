
% constructing a new trellis
v = 11;  % memory element
trellis_1 = poly2trellis(v+1, 4617);
crc_msg_1 = [1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1];
zero_padded_msg_1 = [crc_msg_1, zeros(1, v)];

coded_crc_msg_1 = convenc(zero_padded_msg_1, trellis_1);
cpp_result_1 = [1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1];
areEqual = isequal(coded_crc_msg_1, cpp_result_1);


v2 = 12;
trellis_2 = poly2trellis(v2+1, 17453);
crc_msg_2 = [1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1];
zero_padded_msg_2 = [crc_msg_2, zeros(1, v2)];

coded_crc_msg_2 = convenc(zero_padded_msg_2, trellis_2);
cpp_result_2 = [1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1];
areEqual_2 = isequal(coded_crc_msg_2, cpp_result_2);