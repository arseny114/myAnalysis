#!/bin/bash

Process=E240_qqHinvi
RecoPATH=/cefs/higgs/liugeliang/CEPC/202501/Production/Hinvi/E240_qqHinvi/Reco

# Склеиваем все ana_*.root в один файл
hadd -f ${RecoPATH}/merged_analysis.root ${RecoPATH}/ana_${Process}_*.root

echo "Merged file: ${RecoPATH}/merged_analysis.root"
