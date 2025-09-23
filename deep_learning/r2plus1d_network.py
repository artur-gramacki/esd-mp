import torch
import torch.nn as nn


class Conv2Plus1D(nn.Module):
    def __init__(
            self,
            in_channels: int = 3,
            out_channels: int = 3,
            kernel_size: tuple = (3, 3, 3),
            stride: tuple = (1, 1, 1),
            padding: tuple = (0, 0, 0),
            dilation: tuple = (1, 1, 1),
            bias: bool = True
    ):
        """
        Pytorch implementation of (2+1)D spatiotemporal convolution block.

        This model decomposes a 3D convolution into two separate operations:
        - A 2D spatial convolution (applied per frame) to capture spatial features.
        - A 1D temporal convolution (applied across frames) to model temporal relationships.

        Args:
            in_channels (int): Number of input channels.
            out_channels (int): Number of output channels.
            kernel_size (tuple): (Temporal, Spatial, Spatial) kernel sizes.
            stride (tuple): (Length, Height, Width) strides.
            padding (tuple): (Length, Height, Width) paddings.
            dilation (tuple): (Length, Height, Width) dilation rates.
            bias (bool): Whether to include a bias term in convolutions.

        Attributes:
            spatial_convolution (nn.Conv2d): 2D convolution for spatial feature extraction.
            temporal_convolution (nn.Conv1d): 1D convolution for temporal feature extraction.
            relu (nn.ReLU): Activation function.
        """
        super().__init__()
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.hidden_channels = self.calculate_hidden_channels(kernel_size)
        
        self.spatial_convolution = nn.Conv2d(
            in_channels=self.in_channels,
            out_channels=self.hidden_channels,
            kernel_size=kernel_size[1:],
            stride=stride[1:],
            dilation=dilation[1:],
            padding=padding[1:],
            bias=bias)

        self.temporal_convolution = nn.Conv1d(
            in_channels=self.hidden_channels,
            out_channels=out_channels,
            kernel_size=kernel_size[0],
            stride=stride[0],
            padding=padding[0],
            dilation=dilation[0],
            bias=bias)

        self.relu = nn.ReLU(inplace=True)

    def forward(self, x):
        batch, channels, frames, height, width = x.size()
        x = x.permute(0, 2, 1, 3, 4).contiguous()
        x = x.view(batch * frames, channels, height, width)
        x = self.relu(self.spatial_convolution(x))

        _, channels, height, width = x.size()
        x = x.view(batch, frames, channels, height, width)
        x = x.permute(0, 3, 4, 2, 1).contiguous()
        x = x.view(batch * height * width, channels, frames)
        x = self.temporal_convolution(x)

        channels, frames = x.size(1), x.size(2)
        x = x.view(batch, height, width, channels, frames)
        x = x.permute(0, 3, 4, 1, 2).contiguous()

        return x


    def calculate_hidden_channels(self, kernel_size):
        """
        Computes the number of intermediate channels (projection dimension) 
        used between spatial and temporal convolutions in the (2+1)D block.

        The goal is to approximate the parameter count of a full 3D convolution 
        while maintaining computational efficiency.

        Args:
            kernel_size (tuple): (Temporal, Spatial, Spatial) kernel sizes.

        Returns:
            int: Number of intermediate channels.
        """
        temporal = kernel_size[0]
        spatial = (kernel_size[1] + kernel_size[2]) // 2
                                                                                
        return int(
            (temporal * spatial ** 2 * self.in_channels * self.out_channels) / (
                spatial ** 2 * self.in_channels + temporal * self.out_channels
            )
        )
    

class ConvolutionalBlock(nn.Module):
    def __init__(self, in_channels, out_channels, dropout):
        super().__init__()

        self.convolution = nn.Sequential(
            Conv2Plus1D(in_channels, out_channels, kernel_size=(3, 3, 3), stride=(2, 2, 2), padding=(1, 1, 1)),
            nn.BatchNorm3d(out_channels),
            nn.ReLU(inplace=True),

            Conv2Plus1D(out_channels, out_channels, kernel_size=(3, 3, 3), stride=(1, 1, 1), padding=(1, 1, 1)),
            nn.BatchNorm3d(out_channels),
            nn.Dropout(dropout))
        
        self.shortcut = Conv2Plus1D(in_channels, out_channels, kernel_size=(1, 1, 1), stride=(2, 2, 2), padding=(0, 0, 0))

        self.relu = nn.ReLU(inplace=True)


    def forward(self, x):
        return self.relu(self.convolution(x) + self.shortcut(x))
    

class IdentityBlock(nn.Module):
    def __init__(self, channels, dropout):
        super().__init__()

        self.convolution = nn.Sequential(
            Conv2Plus1D(channels, channels, kernel_size=(3, 3, 3), stride=(1, 1, 1), padding=(1, 1, 1)),
            nn.BatchNorm3d(channels),
            nn.ReLU(inplace=True),

            Conv2Plus1D(channels, channels, kernel_size=(3, 3, 3), stride=(1, 1, 1), padding=(1, 1, 1)),
            nn.BatchNorm3d(channels),
            nn.Dropout(dropout))
        
        self.relu = nn.ReLU(inplace=True)


    def forward(self, x):
        return self.relu(self.convolution(x) + x)
    

class ResidualBlock(nn.Module):
    def __init__(self, in_channels, out_channels, dropout):
        super().__init__()

        self.conv_block = ConvolutionalBlock(in_channels, out_channels, dropout)
        self.identity_block = IdentityBlock(out_channels, dropout)


    def forward(self, x):
        return self.identity_block(self.conv_block(x))


class R2Plus1DConvNet(nn.Module):
    def __init__(
            self,
            in_channels: int = 1,
            num_classes: int = 2,
            dropout: float = 0
    ):
        """
        Pytorch implementation of the R(2+1)D convolutional model for spatiotemporal data
        based on A Closer Look at Spatiotemporal Convolutions for Action Recognition

        This model is designed for video classification and utilizes (2+1)D
        convolutions, which decompose 3D convolutions into separate spatial and
        temporal operations.  The network consists of multiple convolutional blocks
        followed by adaptive pooling  and a fully connected layer for classification.

        Args:
            in_channels (int): Number of channels in the input image.
            num_class (int): Number of output classes for classification.
            dropout (float): : The probability of dropping units during training.
                                                                                
        Attributes:
            net (nn.Sequential): The main body of the network, consisting of:
                - Multiple Conv2Plus1D blocks with LeakyReLU activations.
                - Adaptive average pooling to reduce feature dimensions.
                - Flattening operation to convert tensors to vectors.
                - Fully connected (fc) layer for classification.
        """
        super().__init__()

        self.conv1 = nn.Sequential(
            Conv2Plus1D(in_channels, 16, (3, 7, 7), (1, 2, 2), (1, 3, 3)),
            nn.BatchNorm3d(16),
            nn.ReLU(inplace=True))
        
        self.conv2 = ResidualBlock(16, 16, dropout)
        self.conv3 = ResidualBlock(16, 32, dropout)
        self.conv4 = ResidualBlock(32, 64, dropout)

        self.avgpool = nn.AdaptiveAvgPool3d((1, 1, 1))

        self.fc = nn.Sequential(
            nn.Linear(64, 10),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout),
            nn.Linear(10, num_classes))


    def forward(self, x):
        x = self.conv1(x)
        x = self.conv2(x)
        x = self.conv3(x)
        x = self.conv4(x)

        x = self.avgpool(x)
        x = torch.flatten(x, start_dim=1)

        return self.fc(x)