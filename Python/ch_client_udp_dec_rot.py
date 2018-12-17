import socket
import sys
import math
import struct
import time

numTxAntennas1 = 2
numTxAntennas2 = 0

# Create the datagram socket - Multicast
# sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
# server_address1 = ('192.168.31.12', 50001)
server_address1 = ('10.75.11.164', 50001)
# server_address2 = ('192.168.31.12', 50002)
server_address2 = ('10.75.11.164', 50002)
# server_address2 = ('192.168.31.49', 50002)
# sock.settimeout(20.0)
# print >>sys.stderr, 'connecting to %s port %s' % server_address

# Create the datagram socket - Multicast
sock1 = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
sock2 = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

# Multicast
# multicast_group = ('224.3.29.71', 50001)
# sock.settimeout(20)
# ttl = struct.pack('b', 3)
# sock.setsockopt(socket.IPPROTO_IP, socket.IP_MULTICAST_TTL, ttl)

buffMax1 = 8*2*numTxAntennas1  # Define packet segmentation in characters (characters in TX/RX strings)
buffMax2 = 8*2*numTxAntennas2  # Define packet segmentation in characters (characters in TX/RX strings)
buffMax3 = 8  # One double to synchronize transmissions in time
buffMax4 = 8*2*numTxAntennas1  # One double to send the rotation value to the antennas
buffMax5 = 8*2*numTxAntennas2  # One double to send the rotation value to the antennas
fileName = "helperMUBeamformfeedback1.bin"

data_old1 = ""
data_old2 = ""
counter1 = 0
counter2 = 0

while True:
    f = open(fileName, "rb");

    # Read data
    data1 = f.read(buffMax1)  # 1st transmitter, 2 antennas, 2*8 Bytes each
    #print >> sys.stderr, "1st: %s" % data1.encode('hex_codec')
    data12 = f.read(buffMax3)  # 1st transmitter, 2 antennas, 8 Bytes
    #print >> sys.stderr, "1st: %s" % data12.encode('hex_codec')
    data13 = f.read(buffMax4)  # 1st transmitter, 2 antennas, 8 Bytes
    #print >> sys.stderr, "1st: %s" % data13.encode('hex_codec')
    data2 = f.read(buffMax2)  # 2nd transmitter, 2 antennas, 2*8 Bytes each
    #print >> sys.stderr, "2nd(1): %s" % data2.encode('hex_codec')
    data22 = f.read(buffMax3)  # 2nd transmitter, 1 double, 8 Bytes
    #print >> sys.stderr, "2nd(2): %s" % data22.encode('hex_codec')
    data23 = f.read(buffMax5)  # 1st transmitter, 2 antennas, 8 Bytes
    #print >> sys.stderr, "1st: %s" % data23.encode('hex_codec')

    if data1 == data_old1:
	counter1 = counter1 + 1
    else:
	counter1 = 0

    if data1!=data_old1 or counter1 < 5:
	# Send data
	print >> sys.stderr, '1st sending: "%s"' % data1.encode('hex_codec')
	sock1.sendto(data1+data12+data13, server_address1)
	# Assign new data
	data_old1 = data1

    time.sleep(0.001)

    if data2 == data_old2:
	counter2 = counter2 + 1
    else:
	counter2 = 0

    if data2 and counter2 < 5:
	# Send data
	print >> sys.stderr, '2nd sending: "%s"' % data2.encode('hex_codec')
	time.sleep(0.001)
	sock2.sendto(data2+data22+data23, server_address2)
	# Assign new data
	data_old2 = data2

    f.close()

print >> sys.stderr, 'closing socket'
sock1.close()
sock2.close()
