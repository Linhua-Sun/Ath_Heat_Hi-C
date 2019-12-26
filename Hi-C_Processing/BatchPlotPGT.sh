mkdir  34KEEQC_$(date +%F-%H-%M)
LOG="34KEEQC_$(date +%F-%H-%M).log"
nohup pgt --width 20 --tracks TEKEE34_MinMax.ini --region Chr1:7010662-7110662    -o ./34KEEQC_$(date +%F-%H-%M)/KEE1.pdf  >   ${LOG} 2>&1 &
nohup pgt --width 20 --tracks TEKEE34_MinMax.ini --region Chr2:4105778-4205778    -o ./34KEEQC_$(date +%F-%H-%M)/KEE2.pdf  >>  ${LOG} 2>&1 &
nohup pgt --width 20 --tracks TEKEE34_MinMax.ini --region Chr3:1910791-2010791    -o ./34KEEQC_$(date +%F-%H-%M)/KEE3.pdf  >>  ${LOG} 2>&1 &
nohup pgt --width 20 --tracks TEKEE34_MinMax.ini --region Chr3:3060728-3160728    -o ./34KEEQC_$(date +%F-%H-%M)/KEE4.pdf  >>  ${LOG} 2>&1 &
nohup pgt --width 20 --tracks TEKEE34_MinMax.ini --region Chr3:16652448-16752448  -o ./34KEEQC_$(date +%F-%H-%M)/KEE5.pdf  >>  ${LOG} 2>&1 &
nohup pgt --width 20 --tracks TEKEE34_MinMax.ini --region Chr3:22502744-22602744  -o ./34KEEQC_$(date +%F-%H-%M)/KEE6.pdf  >>  ${LOG} 2>&1 &
nohup pgt --width 20 --tracks TEKEE34_MinMax.ini --region Chr4:11036128-11207665  -o ./34KEEQC_$(date +%F-%H-%M)/KEE7.pdf  >>  ${LOG} 2>&1 &
nohup pgt --width 20 --tracks TEKEE34_MinMax.ini --region Chr4:15439483-15539483  -o ./34KEEQC_$(date +%F-%H-%M)/KEE8.pdf  >>  ${LOG} 2>&1 &
nohup pgt --width 20 --tracks TEKEE34_MinMax.ini --region Chr5:4731440-4831440    -o ./34KEEQC_$(date +%F-%H-%M)/KEE9.pdf  >>  ${LOG} 2>&1 &
nohup pgt --width 20 --tracks TEKEE34_MinMax.ini --region Chr5:10272113-10372113  -o ./34KEEQC_$(date +%F-%H-%M)/KEE10.pdf >>  ${LOG} 2>&1 &




mkdir  WithRepsKEEQC_$(date +%F-%H-%M)
LOG="WithRepsKEEQC_$(date +%F-%H-%M).log"
nohup pgt --width 20 --tracks RepHIC.ini --region Chr1:7010662-7110662    -o ./WithRepsKEEQC_$(date +%F-%H-%M)/KEE1.pdf  >   ${LOG} 2>&1 &
nohup pgt --width 20 --tracks RepHIC.ini --region Chr2:4105778-4205778    -o ./WithRepsKEEQC_$(date +%F-%H-%M)/KEE2.pdf  >>  ${LOG} 2>&1 &
nohup pgt --width 20 --tracks RepHIC.ini --region Chr3:1910791-2010791    -o ./WithRepsKEEQC_$(date +%F-%H-%M)/KEE3.pdf  >>  ${LOG} 2>&1 &
nohup pgt --width 20 --tracks RepHIC.ini --region Chr3:3060728-3160728    -o ./WithRepsKEEQC_$(date +%F-%H-%M)/KEE4.pdf  >>  ${LOG} 2>&1 &
nohup pgt --width 20 --tracks RepHIC.ini --region Chr3:16652448-16752448  -o ./WithRepsKEEQC_$(date +%F-%H-%M)/KEE5.pdf  >>  ${LOG} 2>&1 &
nohup pgt --width 20 --tracks RepHIC.ini --region Chr3:22502744-22602744  -o ./WithRepsKEEQC_$(date +%F-%H-%M)/KEE6.pdf  >>  ${LOG} 2>&1 &
nohup pgt --width 20 --tracks RepHIC.ini --region Chr4:11036128-11207665  -o ./WithRepsKEEQC_$(date +%F-%H-%M)/KEE7.pdf  >>  ${LOG} 2>&1 &
nohup pgt --width 20 --tracks RepHIC.ini --region Chr4:15439483-15539483  -o ./WithRepsKEEQC_$(date +%F-%H-%M)/KEE8.pdf  >>  ${LOG} 2>&1 &
nohup pgt --width 20 --tracks RepHIC.ini --region Chr5:4731440-4831440    -o ./WithRepsKEEQC_$(date +%F-%H-%M)/KEE9.pdf  >>  ${LOG} 2>&1 &
nohup pgt --width 20 --tracks RepHIC.ini --region Chr5:10272113-10372113  -o ./WithRepsKEEQC_$(date +%F-%H-%M)/KEE10.pdf >>  ${LOG} 2>&1 &
