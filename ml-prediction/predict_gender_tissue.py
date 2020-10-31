import numpy as np
import pandas as pd
from model import MLP
from utils import *
import torch
import os, sys, argparse
from Gene import Gene
import datetime
import torch.nn.functional as F
from torch.autograd import Variable
from sklearn.metrics import confusion_matrix, f1_score, accuracy_score
torch.multiprocessing.set_sharing_strategy('file_system')
import seaborn as sns
from sklearn.manifold import TSNE
import pickle as pkl

parser = argparse.ArgumentParser()
parser.add_argument('--n_epoch', type=int, default=50)
parser.add_argument('--seed', type=int, default=1)
parser.add_argument('--print_freq', type=int, default=200)
parser.add_argument('--num_workers', type=int, default=4, help='how many subprocesses to use for data loading')
parser.add_argument('--epoch_decay_start', type=int, default=40)
parser.add_argument('--fold', type=str, default='0')
parser.add_argument('--dataset', type=str, default='GPL11154')
parser.add_argument('--dropout', type=float, default=0.2)
args = parser.parse_args()

# Seed
def seed_everything(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)

    os.environ['PYTHONHASHSEED'] = str(seed)
    torch.backends.cudnn.deterministic = True


# path
dataset = args.dataset
fold = args.fold
dropout_rate = args.dropout
patience = 10

# Hyper Parameters
batch_size = 128
learning_rate = 0.001
if dataset=='GPL570':
    input_size = 135
    out_dim1 = 2
    out_dim2 = 14
elif dataset in ['GPL16791', 'GPL11154', 'GPL18573']:
    input_size = 138
    out_dim1 = 2
    out_dim2 = 14
elif dataset == 'OCTAD':
    input_size = 136
    out_dim1 = 2
    out_dim2 = 13


# load dataset
train_dataset = Gene(root='./input/', dataset=dataset, fold=fold, train=True)
num_iter_per_epoch = int(len(train_dataset.train_labels1)/batch_size)
print('num_iter_per_epoch {:d}'.format(num_iter_per_epoch))
test_dataset = Gene(root='./input/', dataset=dataset, fold=fold, train=False)


# Adjust learning rate and betas for Adam Optimizer
mom1 = 0.9
mom2 = 0.1
alpha_plan = [learning_rate] * args.n_epoch
beta1_plan = [mom1] * args.n_epoch
for i in range(args.epoch_decay_start, args.n_epoch):
    alpha_plan[i] = float(args.n_epoch - i) / (args.n_epoch - args.epoch_decay_start) * learning_rate
    beta1_plan[i] = mom2

def adjust_learning_rate(optimizer, epoch):
    for param_group in optimizer.param_groups:
        param_group['lr']=alpha_plan[epoch]
        param_group['betas']=(beta1_plan[epoch], 0.999) # Only change beta1


save_dir = './output/' + dataset + '_task_2/'

if not os.path.exists(save_dir):
    os.system('mkdir -p %s' % save_dir)

txtfile = save_dir + '/dropout_' + str(dropout_rate) + '_fold_' + str(fold) + '.txt'
nowTime = datetime.datetime.now().strftime('%Y-%m-%d-%H:%M:%S')
if os.path.exists(txtfile):
    os.system('mv %s %s' % (txtfile, txtfile+".bak-%s" % nowTime))  # rename exsist file for collison


def masked_accuracy(logit, target):
    unmask_idx = np.where(target.cpu() != -1)[0]

    output = F.softmax(logit, dim=1)
    _, pred = torch.max(output.data, 1)

    acc = (pred[unmask_idx].cpu() == target[unmask_idx].cpu()).sum() / float(len(unmask_idx))

    return acc


def loss_fun(logits, labels):
    mask = labels.cpu()==-1
    loss = F.cross_entropy(logits[~mask], labels[~mask].long())
    return loss


def train(train_loader, epoch, model, optimizer):
    train_total1 = 0
    train_correct1 = 0
    train_total2 = 0
    train_correct2 = 0

    for i, (images, labels1, labels2, indexes) in enumerate(train_loader):
        if i > num_iter_per_epoch:
            break

        images = Variable(images).cuda()
        labels1 = Variable(labels1).cuda()
        labels2 = Variable(labels2).cuda()

        # Forward + Backward + Optimize
        logits1, logits2, _, _ = model(images)
        acc1 = masked_accuracy(logits1, labels1)
        acc2 = masked_accuracy(logits2, labels2)
        train_total1 += 1
        train_correct1 += acc1
        train_total2 += 1
        train_correct2 += acc2

        loss1 = loss_fun(logits1, labels1)
        loss2 = loss_fun(logits2, labels2)
        loss = loss1 + loss2

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        if (i + 1) % args.print_freq == 0:
            print('epoch {:d}/{:d} Iter {:d}/{:d} trn_acc1 {:.4f} trn_acc2 {:.4f}'.format(epoch, args.n_epoch, i+1, num_iter_per_epoch, acc1, acc2))

    train_acc1 = float(train_correct1) / float(train_total1)
    train_acc2 = float(train_correct2) / float(train_total2)
    return train_acc1, train_acc2


def eval_metrics(labels_list, preds_list):
    """list of batch labels and batch preds"""
    # faltten first
    labels_flatten = np.array([item.data for sublist in labels_list for item in sublist], dtype=np.int32)
    preds_flatten = np.array([item.data for sublist in preds_list for item in sublist], dtype=np.int32)

    mask = labels_flatten==-1
    #print(mask)

    cm = confusion_matrix(labels_flatten[~mask], preds_flatten[~mask])
    f1 = f1_score(labels_flatten[~mask], preds_flatten[~mask], average='macro')
    acc = accuracy_score(labels_flatten[~mask], preds_flatten[~mask])
    return cm, f1, acc


# Evaluate the Model
def evaluate(test_loader, model, epoch):
    model.eval()    # Change model to 'eval' mode.
    pred_list1 = []
    labels_list1 = []

    pred_list2 = []
    labels_list2 = []
    for images, labels1, labels2, _ in test_loader:
        images = Variable(images).cuda()
        logits1, logits2, _, _ = model(images)

        outputs1 = F.softmax(logits1, dim=1)
        _, pred1 = torch.max(outputs1.data, 1)

        pred_list1.append(pred1.cpu())
        labels_list1.append(labels1)

        outputs2 = F.softmax(logits2, dim=1)
        _, pred2 = torch.max(outputs2.data, 1)

        pred_list2.append(pred2.cpu())
        labels_list2.append(labels2)

    cm1, f11, acc1 = eval_metrics(labels_list1, pred_list1)
    cm2, f12, acc2 = eval_metrics(labels_list2, pred_list2)

    if epoch==49:
        print('acc_1 {:.2f} f1_score {:.2f}'.format(acc1, f11))
        print('acc_2 {:.2f} f1_score {:.2f}'.format(acc2, f12))
        print(cm1)
        print(cm2)

    # cross table for gender_pred and tissue_pred
    print(pd.crosstab(np.array([item.data for sublist in pred_list1 for item in sublist], dtype=np.int32), np.array([item.data for sublist in pred_list2 for item in sublist], dtype=np.int32)))
    return acc1, f11, acc2, f12


def get_representation(test_loader, model):
    model.eval()    # Change model to 'eval' mode.
    rep1_matrix = np.empty([0, 32], dtype=np.float32)
    rep2_matrix = np.empty([0, 32], dtype=np.float32)
    label1_list = []
    label2_list = []
    for images, labels1, labels2, _ in test_loader:
        images = Variable(images).cuda()
        _, _, rep1, rep2 = model(images)
        rep1_matrix = np.concatenate([rep1_matrix, rep1.cpu().detach().numpy()])
        rep2_matrix = np.concatenate([rep2_matrix, rep2.cpu().detach().numpy()])

        label1_list += labels1
        label2_list += labels2

    return rep1_matrix, rep2_matrix, label1_list, label2_list


def output_prediction(test_loader, model):
    model.eval()    # Change model to 'eval' mode.
    OUT1 = np.empty([0, out_dim1], dtype=np.float32)
    OUT2 = np.empty([0, out_dim2], dtype=np.float32)
    for images, labels1, labels2, _ in test_loader:
        images = Variable(images).cuda()
        logits1, logits2, _, _ = model(images)

        outputs1 = F.softmax(logits1, dim=1).cpu().detach().numpy()
        outputs2 = F.softmax(logits2, dim=1).cpu().detach().numpy()

        OUT1 = np.concatenate([OUT1, outputs1], axis=0)
        OUT2 = np.concatenate([OUT2, outputs2], axis=0)

    OUT = np.concatenate([OUT1, OUT2], axis=1)
    fn = save_dir + 'fold_all_predictions.txt'
    np.savetxt(fn, OUT, delimiter='\t')
    print(OUT.shape)
    print(test_dataset.test_data.shape)


def main():
    seed_everything(args.seed)
    # Data Loader (Input Pipeline)
    train_loader = torch.utils.data.DataLoader(dataset=train_dataset,
                                               batch_size=batch_size,
                                               num_workers=args.num_workers,
                                               drop_last=True,
                                               shuffle=True)

    test_loader = torch.utils.data.DataLoader(dataset=test_dataset,
                                              batch_size=batch_size,
                                              num_workers=args.num_workers,
                                              drop_last=False,
                                              shuffle=False)

    # Define models
    mlp = MLP(input_size=input_size, out_dim1=out_dim1, out_dim2=out_dim2, dropout_rate=dropout_rate)
    mlp.cuda()
    print(mlp.parameters)
    optimizer1 = torch.optim.Adam(mlp.parameters(), lr=learning_rate)

    with open(txtfile, "a") as myfile:
        myfile.write('epoch: train_acc1 train_acc2 test_acc1 test_acc2 test_f1 test_f2 \n')

    epoch = 0
    train_acc1 = 0
    train_acc2 = 0
    # evaluate models with random weights
    test_acc1, test_f11, test_acc2, test_f12 = evaluate(test_loader, mlp, epoch)
    print('Epoch {:d}/{:d} tst_acc1 {:.4f} tst_acc2 {:.4f} tst_f11 {:.4f} tst_f12 {:.4f}'.format(epoch, args.n_epoch, test_acc1, test_acc2, test_f11, test_f12))
    # save results
    with open(txtfile, "a") as myfile:
        myfile.write(str(int(epoch)) + ': ' + str(train_acc1) + ' ' + str(train_acc2) + ' ' + str(test_acc1) + ' '
                     + str(test_acc2) + ' ' + str(test_f11) + ' ' +  str(test_f12) + ' '  + "\n")

    # training
    best_acc = 0.0
    wait = 0
    for epoch in range(1, args.n_epoch):
        # train models
        mlp.train()
        adjust_learning_rate(optimizer1, epoch)
        train_acc1, train_acc2 = train(train_loader, epoch, mlp, optimizer1)
        # evaluate models
        test_acc1, test_f11, test_acc2, test_f12 = evaluate(test_loader, mlp, epoch)
        # save results
        print(
            'Epoch {:d}/{:d} tst_acc1 {:.4f} tst_acc2 {:.4f} tst_f11 {:.4f} tst_f12 {:.4f}'.format(epoch, args.n_epoch,
                                                                                                   test_acc1, test_acc2,
                                                                                                   test_f11, test_f12))
        with open(txtfile, "a") as myfile:
            myfile.write(str(int(epoch)) + ': ' + str(train_acc1) + ' ' + str(train_acc2) + ' ' + str(test_acc1) + ' '
                         + str(test_acc2) + ' ' + str(test_f11) + ' ' + str(test_f12) + ' ' + "\n")

        # Early stopping
        if np.greater(test_acc1+test_acc2, best_acc):
            best_acc = test_acc1+test_acc2
            wait = 0
        else:
            wait += 1
            if wait >= patience:
                print('Early stop at Epoch: {:d} with final val acc: {:.4f} {:.4f}'.format(epoch, test_acc1, test_acc2))
                rep1_matrix, rep2_matrix, label1_list, label2_list = get_representation(test_loader, mlp)
                fn = save_dir + 'rep_' + fold + '.pkl'
                with open(fn, 'wb') as f:
                    pkl.dump([rep1_matrix, rep2_matrix, label1_list, label2_list], f)
                if fold == 'all':
                    output_prediction(test_loader, mlp)
                break

    if fold == 'all':
        output_prediction(test_loader, mlp)


if __name__ == '__main__':
    main()
