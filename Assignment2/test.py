import torch
dtype = torch.float
device = torch.device("mps") # This executes all calculations on the CPU
# device = torch.device("cuda:0") # This executes all calculations on the GPU

# Creation of a tensor and filling of a tensor with random numbers
a = torch.randn(2, 3, device=device, dtype=dtype)
print(a) # Output of tensor A
# Output: tensor([[-1.1884,  0.8498, -1.7129],
#                  [-0.8816,  0.1944,  0.5847]])

# Creation of a tensor and filling of a tensor with random numbers
b = torch.randn(2, 3, device=device, dtype=dtype)
print(b) # Output of tensor B
# Output: tensor([[ 0.7178, -0.8453, -1.3403],
#                  [ 1.3262,  1.1512, -1.7070]])

print(a*b) # Output of a multiplication of the two tensors
# Output: tensor([[-0.8530, -0.7183,  2.58],
#                  [-1.1692,  0.2238, -0.9981]])

print(a.sum()) # Output of the sum of all elements in tensor A
# Output: tensor(-2.1540)

print(a[1,2]) # Output of the element in the third column of the second row (zero based)
# Output: tensor(0.5847)

print(a.max()) # Output of the maximum value in tensor A
# Output: tensor(-1.7129)