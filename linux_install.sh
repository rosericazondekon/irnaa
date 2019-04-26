#!/bin/bash
echo $'export PATH=$PATH:'$(pwd)/linux >> ~/.bashrc && echo $'export IRNAAPATH='$(pwd)/linux  >> ~/.bashrc && source  ~/.bashrc
