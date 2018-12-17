import matlab.engine
import sys
import time
import wx
import os
import ntplib
from gnuradio import uhd
from gnuradio import gr
from PyQt4 import Qt
import subprocess

start1 = time.time()

class top_block(gr.top_block, Qt.QWidget):

    def __init__(self,IPAddress):

        ##################################################
        # Blocks
        ##################################################
        self.uhd_usrp_sink_0 = uhd.usrp_sink(
        	" ,".join((IPAddress, "")),
        	uhd.stream_args(
        		cpu_format="fc32",
        		channels=range(1),
        	),
        )
        self.uhd_usrp_sink_0.set_clock_source("external", 0)
        self.uhd_usrp_sink_0.set_samp_rate(32000)
        self.uhd_usrp_sink_0.set_center_freq(915e6, 0)
        self.uhd_usrp_sink_0.set_gain(0, 0)
        self.uhd_usrp_sink_0.set_antenna("TX/RX", 0)
        #self.analog_sig_source_x_0 = analog.sig_source_c(samp_rate, analog.GR_COS_WAVE, 1000, 1, 0)


	##################################################
	# Adding line for time sync purpose
	##################################################

	try:
	    client = ntplib.NTPClient()
	    response = client.request('pool.ntp.org')
	    os.system('date ' + time.strftime('%m%d%H%M%Y.%S',time.localtime(response.tx_time)))
	    #os.system(time.strftime('%m%d%H%M%Y.%S',time.localtime(response.tx_time)))
	    system_time_min = int(time.strftime('%M'))
	    system_time_min = (system_time_min*60)
            #print (system_time_min)
            system_time_sec = int(time.strftime('%S'))
            #print (system_time_sec)
            system_time_str = (system_time_min + system_time_sec)
	    #print(system_time_str)
            system_time = int(system_time_str)
            print (system_time)
	except:
	    print('Could not sync with time server.')

	print('Done.')
	#time.sleep(2.0)
	#self.uhd_usrp_source_0.set_time_unknown_pps(uhd.time_spec(0.0))
        #self.uhd_usrp_source_0.set_time_next_pps(system_time)
	#self.uhd_usrp_source_0.set_time_now(uhd.time_spec(system_time))
        self.uhd_usrp_sink_0.set_time_next_pps(uhd.time_spec_t(system_time +1))
	##################################################

	end1 = time.time()

	print "radio instantiation time is ", end1 - start1, " seconds"

        ##################################################
        # Getting usrp time
        ##################################################
        for ii in range(0,5):
        	#last_time = self.uhd_usrp_source_0.get_time_now()
		last_pps0 = self.uhd_usrp_sink_0.get_time_last_pps()
		#print "last_time : %6.12f\n"%uhd.time_spec_t.get_real_secs(last_time)
		print "last_time_0 : %6.12f\n"%uhd.time_spec_t.get_real_secs(last_pps0)
		#print "last_pps0 : %6.12f\n"%uhd.time_spec_t.get_real_secs(last_pps0)
		time.sleep(1.0)
	##################################################


def waitForOthers(t):
	################################################## 
	# Wwaiting for all radios to sync in time
	##################################################
	print "waiting for other radios to get synced\n"

	################################################## 
	# Waiting for all radios to sync in time
	##################################################
	check_time = time.time()
	warn = check_time - start1

	print "time elapsed: ", warn

	if warn > t:
		print "time elapsed > wait time given. system will go into infinite loop: ", warn

	while True:
		end_time = time.time()
		m = int(end_time - start1)
		print m
		if m == t:
			break;

def initMatlabEngine():
	################################################## 
	# Initializing MATLAB engine
	##################################################
	print "instantiating matlab engine\n"
	myMatlab = matlab.engine.start_matlab()	

	return myMatlab



def main(top_block_cls=top_block, options=None):

#	IPAddress = "addr=192.168.10.2"
#	top_block_cls(IPAddress)

	myMatlab1 = initMatlabEngine()

	waitForOthers(30)

	myMatlab1.BER_MultiUserBeamformingExample4_tx1(nargout=0)



if __name__ == '__main__':
    main()
