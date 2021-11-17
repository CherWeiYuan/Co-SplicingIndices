# Co-Splicing Indices (Prototype, requires stress test)

Co-splicing Indices gives a quantitative measure of the strength and confidence of co-splicing.

Invitees to this repository are given a detailed powerpoint to explain the the mathematical and algorithmic concept of CSI.

## Installations
PRIMERg runs on Linux. For Windows users, you can download Ubuntu (tested on Ubutun 20.04 but new versions should work. Installation guide: https://ubuntu.com/tutorials/ubuntu-on-windows#1-overview)

First, run: ```sudo apt update``` and ```sudo apt upgrade```

Check if you have Python 3.8.10 on Ubuntu (Confirm that Python3 is already install with ```python3 --version```. If you need to update your python version, do: ```sudo apt update```, ```sudo apt upgrade```, then ```sudo apt upgrade python3```)

Get pip for python: ```sudo apt install python3-pip```

Get required python packages:
```pip install plotly==5.3.1```
```pip install scipy==1.7.2```
```pip install pandas==1.3.4```

## Quickstart
You can try CSI.py with sample test_input.csv [here](https://github.com/CherWeiYuan/Co-Splicing_Indices/tree/main/sample_input_output)
1. Create a new folder. Put CSI.py and user_inputs.py in new folder.
2. Edit user_inputs.py using any text editor.
3. Navigate to new folder in Ubuntu. For example if your new folder named 'CSI' is on Desktop ```cd /mnt/c/Users/{username}/Desktop/CSI```
4. Run ```python3 csi.py``` on Ubuntu.

