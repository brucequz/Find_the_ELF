
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


% Create two row vectors
vector1 = [1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1];
vector2 = [1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1];

% Calculate the total length of the interleaved vector
totalLength = numel(vector1) + numel(vector2);

% Initialize the interleaved vector with zeros
interleavedVector = zeros(1, totalLength);

% Interleave the elements from both vectors
interleavedVector(2:2:end) = vector1;
interleavedVector(1:2:end) = vector2;

% Display the interleaved vector
disp(interleavedVector);

cpp_result_interleave = [1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1];

areEqual_3 = isequal(interleavedVector, cpp_result_interleave);

crc_poly = [1, 1, 0, 1];
convolved = conv([1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1],[1, 1, 0, 1]);
convolved = mod(convolved, 2);

cpp_result_4 = [1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1];
areEqual_4 = isequal(convolved, cpp_result_4);

[Q, R] = deconv(cpp_result_4, crc_poly);
Q = mod(Q, 2);
