INCFLAGS = 
LDFLAGS =
LIBLIST = -lm

CC = gcc
CXX = g++
######## config part ########
ResultsDir="/work/data2/yamada/OmgOmg/Set2/Spin0-0"
#ResultsDir="./hoge"
SrcSpin=0
SinkSpin=0
######## config part ########
#BSUB = bsub
#LOG=  -o ../log/$@
TDER = Tderlog.out
#TDER = Tdernonlog.out

################
#2point
SRC_OmgOmg2pt=OmgOmg2pt.cpp                                                                                                                                                                                                
SRC_jack_sub2point=jack_sub2point.cpp
SRC_jack_2point.tcorr.final=jack_2point.tcorr.final.cpp
SRC_jack_2pointfinal=jack_2pointfinal.cpp
#4point bin & analys                                                                                                                                                                                   
SRC_jack_subBSwave=jack_subBSwave.cpp
SRC_makeR=makeR.cpp
SRC_Laplacian=Laplacian.cpp
SRC_Tderlog=Tderlog.cpp
SRC_Tder=Tdernonlog.cpp
SRC_sum=sum.cpp
#jack final results                                                                                                                                                                                    
SRC_jack_BSfinal=jack_BSfinal.cpp
SRC_jack_Rfinal=jack_Rfinal.cpp
SRC_jack_Tderfinal=jack_Tderfinal.cpp
SRC_jack_Lapfinal=jack_Lapfinal.cpp
SRC_jack_final=jack_final.cpp
SRC_aveTshift=aveTshift.cpp


#2point
TARGET_jack_sub2point=jack_sub2point.out
TARGET_jack_2point.tcorr.final=jack_2point.tcorr.final.out
TARGET_jack_2pointfinal=jack_2pointfinal.out
#4point bin & analys
TARGET_jack_subBSwave=jack_subBSwave.out
TARGET_makeR=makeR.out
TARGET_Laplacian=Laplacian.out
TARGET_Tderlog=Tderlog.out
TARGET_Tder=Tdernonlog.out
TARGET_sum=sum.out
#jack final results
TARGET_jack_BSfinal=jack_BSfinal.out
TARGET_jack_Rfinal=jack_Rfinal.out
TARGET_jack_Tderfinal=jack_Tderfinal.out
TARGET_jack_Lapfinal=jack_Lapfinal.out
TARGET_jack_final=jack_final.out
TARGET_aveTshift=aveTshift.out


TARGET_LIST = $(TARGET_OmgOmg2pt) \
$(TARGET_jack_sub2point) \
$(TARGET_jack_2point.tcorr.final) \
$(TARGET_jack_subBSwave) \
$(TARGET_makeR) \
$(TARGET_Laplacian) \
$(TARGET_Tderlog) \
$(TARGET_Tder) \
$(TARGET_sum) \
$(TARGET_jack_BSfinal) \
$(TARGET_jack_Rfinal) \
$(TARGET_jack_Tderfinal) \
$(TARGET_jack_Lapfinal) \
$(TARGET_jack_final) \
$(TARGET_aveTshift)

################ make part ###########################
all:$(TARGET_LIST)

$(TARGET_jack_sub2point): $(SRC_jack_sub2point)
	$(CXX) $(LDFLAGS) -o $@ $(SRC_jack_sub2point)
$(TARGET_jack_2point.tcorr.final): $(SRC_jack_2point.tcorr.final)	
	$(CXX) $(LDFLAGS) -o $@ $(SRC_jack_2point.tcorr.final)
$(TARGET_jack_2pointfinal): $(SRC_jack_2pointfinal)	
	$(CXX) $(LDFLAGS) -o $@ $(SRC_jack_2pointfinal)
$(TARGET_jack_subBSwave): $(SRC_jack_subBSwave)	
	$(CXX) $(LDFLAGS) -o $@ $(SRC_jack_subBSwave)
$(TARGET_makeR): $(SRC_makeR)	
	$(CXX) $(LDFLAGS) -o $@ $(SRC_makeR)
$(TARGET_Laplacian): $(SRC_Laplacian)	
	$(CXX) $(LDFLAGS) -o $@ $(SRC_Laplacian)
$(TARGET_Tderlog): $(SRC_Tderlog)	
	$(CXX) $(LDFLAGS) -o $@ $(SRC_Tderlog)
$(TARGET_Tder): $(SRC_Tder)	
	$(CXX) $(LDFLAGS) -o $@ $(SRC_Tder)
$(TARGET_sum): $(SRC_sum)	
	$(CXX) $(LDFLAGS) -o $@ $(SRC_sum)
$(TARGET_jack_BSfinal): $(SRC_jack_BSfinal)	
	$(CXX) $(LDFLAGS) -o $@ $(SRC_jack_BSfinal)
$(TARGET_jack_Rfinal): $(SRC_jack_Rfinal)	
	$(CXX) $(LDFLAGS) -o $@ $(SRC_jack_Rfinal)
$(TARGET_jack_Tderfinal): $(SRC_jack_Tderfinal)	
	$(CXX) $(LDFLAGS) -o $@ $(SRC_jack_Tderfinal)
$(TARGET_jack_Lapfinal): $(SRC_jack_Lapfinal)
	$(CXX) $(LDFLAGS) -o $@ $(SRC_jack_Lapfinal)
$(TARGET_jack_final): $(SRC_jack_final)
	$(CXX) $(LDFLAGS) -o $@ $(SRC_jack_final)
$(TARGET_aveTshift): $(SRC_aveTshift)
	$(CXX) $(LDFLAGS) -o $@ $(SRC_aveTshift)


clean:
	rm -rf ./*.out ./*~
logclean:
	rm ../log/*
init:
	mkdir -p $(ResultsDir)/Projwave/xyz/Spin$(SinkSpin)-$(SrcSpin)
	mkdir -p $(ResultsDir)/tcor2pt/raw
	@echo "----------------------------------------------------------"
	@echo "you should copy input data to "
	@echo " projected BSwave: $(ResultsDir)/Projwave/xyz/Spin$(SinkSpin)-$(SrcSpin)"
	@echo " time correlator :  $(ResultsDir)/tcor2pt/raw"
	@echo "----------------------------------------------------------"
#remain ProjBSwave
dataclean:
	rm -r $(ResultsDir)/Laplacian
	rm -r $(ResultsDir)/Potential
	rm -r $(ResultsDir)/impBSwave
	rm -r $(ResultsDir)/Tder
	rm -r $(ResultsDir)/Rcor
	rm -r $(ResultsDir)/Projwave/binProjwave
	rm -r $(ResultsDir)/Projwave/jack_error
	rm -r $(ResultsDir)/tcor2pt/bintcor2pt
	rm -r $(ResultsDir)/tcor2pt/jack_error


help:
	@echo "#####################################################"
	@echo "###                how to use                     ###"
	@echo "#####################################################"
	@echo "------- data stream--------------------------------"
	@echo "noise_redection=>aveTshift=>js4pt"
	@echo "2pt=>js2pt"
	@echo "js4pt,js2pt>binRcor=>binLaplacian,binTder=>Potential,fit"
	@echo "------- comand naming memo--------------------------------"
	@echo "js: jack sub (bin average)	"
	@echo "jf: jack final results (averageg , error)"
	@echo "	"
	@echo "-------- clean command ---------------------"
	@echo "clean    : clean *.out *.~"
	@echo "logclean : clean logfile"
	@echo " "
	@echo "-------- config memo ----------------------"
	@echo " You should change only config part[output path & some swtich] in Makefile"
	@echo " You should change only some parametor in include/analys.h"
	@echo "TDER = Tderlog.out     : -d/dt(logR(t))"
	@echo "TDER = Tdernonlog.out  : -(1/R(t))*d/dt(R(t))"
	@echo "SrcSpin                :Src  Spin"
	@echo "SinkSpin               :Sink Spin"
	@echo " "
	@echo "-------- directory info ----------------------"
	@echo "ResultsDir---Laplacian --binLaplacian---xyz"
	@echo "           |           |-jack_error  ---xyz"
	@echo "           |                          |--r"
	@echo "           |"
	@echo "           |-Potential --binPotential---xyz"
	@echo "           |           |"
	@echo "           |           |-jack_error ---xyz"
	@echo "           |           |               |---r"
	@echo "           |"
	@echo "           |-impBSwave --aveTshift -xyz"
	@echo "           |           |-noize_redection--ts{}"
	@echo "           |                          "
	@echo "           |"
	@echo "           |-Projwave ---binProjwave ---xyz"
	@echo "           |           |-jack_error  ---xyz"
	@echo "           |           |              |--r"
	@echo "           |           |-xyz---Spin"
	@echo "           |"
	@echo "           |-Rcor -------binR ----------xyz"
	@echo "           |           |-jack_error  ---xyz"
	@echo "           |                          |--r"
	@echo "           |"
	@echo "           |-tcor2pt ----row"
	@echo "           |           |-bintcor2pt ----xyz"
	@echo "           |           |-jack_error  ---xyz"
	@echo "           |                          |--r"
	@echo "           |"
	@echo "           |-Tder -------binTder -------xyz"
	@echo "                       |-jack_error  ---xyz"
	@echo "                                      |--r"
	@echo "-------- run command -----------------------"
	@echo "runjs2pt  : 2pt time correlator  ==> 2pt time correlator (time shift averaged)==> bin 2pt time correlator"
	@echo "runtave   : 4pt correlator                       ==> 4pt correlator (time shift averaged)"
	@echo "runjsBS   : 4pt correlator (time shift averaged) ==> bin 4pt correlator"
	@echo "runjf2pt  : bin 2pt correlator                   ==> jack 2pt correlator"
	@echo "runmakeR  : 2pt + 4pt (bin corre)                ==> bin R correlator"
	@echo "runLap    : bin R correlator                     ==> bin Laplacian part"
	@echo "runTder   : bin R correlator                     ==> bin Tdepndent part"
	@echo "runsum    : Lap + Tdep (bin)                     ==> bin Potencial"
	@echo "runjfBS   : bin 4pt correlator                   ==> jack 4pt correlator"
	@echo "runjfR    : bin R correlator                     ==> jack R corelator"
	@echo "runjfTder : bin Tdependent part                  ==> jack Tdependent part"
	@echo "runjfLap  : bin Laplacian part                   ==> jack Laplacian part"
	@echo "runjfPot  : bin Potencial                        ==> jack Potencial"
	@echo "#######################################################"
############### run part ##############
runjs2pt:
		$(BSUB) $(LOG) ./$(TARGET_jack_sub2point) $(ResultsDir)
runjf2pt:
		$(BSUB) $(LOG) ./$(TARGET_jack_2point.tcorr.final) $(ResultsDir)
runjsBS:
		$(BSUB) $(LOG) ./$(TARGET_jack_subBSwave) $(ResultsDir)
runmakeR:
		$(BSUB) $(LOG) ./$(TARGET_makeR) $(ResultsDir)
runLap:
		$(BSUB) $(LOG) ./$(TARGET_Laplacian) $(ResultsDir)
runTder:
		$(BSUB) $(LOG) ./$(TDER) $(ResultsDir)
runsum:
		$(BSUB) $(LOG) ./$(TARGET_sum) $(ResultsDir)
runjfBS:
		$(BSUB) $(LOG) ./$(TARGET_jack_BSfinal) $(ResultsDir)
runjfR:
		$(BSUB) $(LOG) ./$(TARGET_jack_Rfinal) $(ResultsDir)
runjfTder:
		$(BSUB) $(LOG) ./$(TARGET_jack_Tderfinal) $(ResultsDir)
runjfLap:
		$(BSUB) $(LOG) ./$(TARGET_jack_Lapfinal) $(ResultsDir)
runjfPot:
		$(BSUB) $(LOG) ./$(TARGET_jack_final) $(ResultsDir)
runtave:
		$(BSUB) $(LOG) ./$(TARGET_aveTshift) $(ResultsDir)
############# rule ################
.cpp.out:
	$(CXX) $(INCFLAGS) -o $@ $<

