import math
import torch
import torch.nn as nn
import torch.nn.init as init 
import torch.nn.functional as F
import torch.optim as optim

def call_bn(bn, x):
    return bn(x)

class MLP(nn.Module):

    def __init__(self, input_size=135, out_dim1=2, out_dim2=14, dropout_rate=0.2, top_bn=False):
        self.dropout_rate = dropout_rate
        self.top_bn = top_bn
        super(MLP, self).__init__()
        self.fc1 = nn.Linear(input_size, 64)
        self.fc2 = nn.Linear(64, 32)

        self.sub1_fc3 = nn.Linear(32, 32)
        self.sub1_fc4 = nn.Linear(32, out_dim1)

        self.sub2_fc3 = nn.Linear(32, 32)
        self.sub2_fc4 = nn.Linear(32, out_dim2)

    def forward(self, x):
        h=x
        h=self.fc1(h)
        h=F.relu(h)
        h = F.dropout2d(h, p=self.dropout_rate)

        h=self.fc2(h)
        h=F.relu(h)
        h = F.dropout2d(h, p=self.dropout_rate)

        sub1_h = self.sub1_fc3(h)
        sub1_h = F.relu(sub1_h)
        logit1 = self.sub1_fc4(sub1_h)

        sub2_h = self.sub2_fc3(h)
        sub2_h = F.relu(sub2_h)
        logit2 = self.sub2_fc4(sub2_h)

        return logit1, logit2, sub1_h, sub2_h


class MLP3(nn.Module):
    def __init__(self, input_size=135, out_dim1=2, out_dim2=14, out_dim3=3, dropout_rate=0.2, top_bn=False):
        self.dropout_rate = dropout_rate
        self.top_bn = top_bn
        super(MLP3, self).__init__()
        self.fc1 = nn.Linear(input_size, 64)
        self.fc2 = nn.Linear(64, 32)

        self.sub1_fc3 = nn.Linear(32, 32)
        self.sub1_fc4 = nn.Linear(32, out_dim1)

        self.sub2_fc3 = nn.Linear(32, 32)
        self.sub2_fc4 = nn.Linear(32, out_dim2)

        self.sub3_fc3 = nn.Linear(32, 32)
        self.sub3_fc4 = nn.Linear(32, out_dim3)

    def forward(self, x):
        h=x
        h=self.fc1(h)
        h=F.relu(h)
        h = F.dropout2d(h, p=self.dropout_rate)

        h=self.fc2(h)
        h=F.relu(h)
        h = F.dropout2d(h, p=self.dropout_rate)

        sub1_h = self.sub1_fc3(h)
        sub1_h = F.relu(sub1_h)
        logit1 = self.sub1_fc4(sub1_h)

        sub2_h = self.sub2_fc3(h)
        sub2_h = F.relu(sub2_h)
        logit2 = self.sub2_fc4(sub2_h)

        sub3_h = self.sub3_fc3(h)
        sub3_h = F.relu(sub3_h)
        logit3 = self.sub3_fc4(sub3_h)

        return logit1, logit2, logit3