#!/bin/bash

objects=(EleInfo EventInfo GenJetInfo GenPartInfo JetInfo MetInfo MhtInfo MuonInfo PairInfo TauInfo VertexInfo)
webdir="https://raw.githubusercontent.com/acarnes/UfHMuMuCode/mu_e_CMSSW_8_0_X_AWB/UFDiMuonsAnalyzer/"

for i in "${objects[@]}"
do
    wget ${webdir}"interface/"$i".h" --no-check-certificate 
    wget ${webdir}"src/"$i".cc" --no-check-certificate 
done

sed -i "s/UfHMuMuCode\/UFDiMuonsAnalyzer\/interface\///g" *.cc
sed -i "s/import/include/g" *.cc
sed -i '/typedef/d' *.h
