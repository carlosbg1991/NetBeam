import socket
import sys
import math
import struct
import time

numTxAntennas = 4

# Create the datagram socket - Multicast
# sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
server_address = ('192.168.31.41', 50001)
# sock.settimeout(20.0)
# print >>sys.stderr, 'connecting to %s port %s' % server_address

# Create the datagram socket - Multicast
sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

# Multicast
# multicast_group = ('224.3.29.71', 50001)
# sock.settimeout(20)
# ttl = struct.pack('b', 3)
# sock.setsockopt(socket.IPPROTO_IP, socket.IP_MULTICAST_TTL, ttl)

buffMax = 8*2*numTxAntennas  # Define packet segmentation in characters (characters in TX/RX strings)
fileName = "helperMUBeamformfeedback1.bin"

data_old = ""
counter = 0

while True:
    f = open(fileName, "rb");

    # Read data
    data = f.read(buffMax)  # 1st transmitter, 2 antennas, 32 Bytes

    if data == data_old:
	counter = counter + 1
    else:
	data_old = data
	counter = 0

    if data and counter < 10:
	# Send data
	print >> sys.stderr, '1st sending: "%s"' % data.encode('hex_codec')
	sock.sendto(data, server_address)

    time.sleep(0.001)

    f.close()

print >> sys.stderr, 'closing socket'
sock1.close()
