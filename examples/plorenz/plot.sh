#!/bin/bash

gnuplot -e "filename='datas' ; output_file='datas.png'" -p -c plot.cfg
gnuplot -e "filename='data_cov1' ; output_file='data_cov1.png'" -p -c plot.cfg
gnuplot -e "filename='data_cov2' ; output_file='data_cov2.png'" -p -c plot.cfg
gnuplot -e "filename='data_cov3' ; output_file='data_cov3.png'" -p -c plot.cfg
gnuplot -e "filename='data_cov_mc' ; output_file='data_cov_mc.png'" -p -c plot.cfg
gnuplot -e "filename='datas_mc' ; output_file='datas__mc.png'" -p -c plot.cfg

