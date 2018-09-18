import socket
import sys
import time
import struct
import statistics
import math
# import numpy
# import matplotlib.pyplot as pp

# Create UDP socket 
sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
server_address = ('', 50001)
sock.bind(server_address)
# sock.setblocking(0)
print >>sys.stderr, 'starting up on %s port %s' % server_address  # Default print

# Multicast communication
# multicast_group = '224.3.29.71'
# group = socket.inet_aton(multicast_group)
# mreq = struct.pack('4sL', group, socket.INADDR_ANY)
# sock.setsockopt(socket.IPPROTO_IP, socket.IP_ADD_MEMBERSHIP, mreq)

buffMax = 64  # Define packet segmentation (characters in TX/RX strings)
name = "mierda.bin"

results = []
mean = 0
updateIdx = 0;
while True:    
    wfd = open(name, "wb")

    data_old = ""

    try:
        while True:
            var1 = time.time()
            data, address = sock.recvfrom(buffMax)

            if data:
               # if data != data_old:
               #     data_old = data
                wfd.seek(0)
                # wfd.truncate()
                wfd.write(data)
                print "Received data: %s" % data.encode('hex_codec')
                var2 = time.time()
                updateIdx = updateIdx + 1;
                # results.append(var2-var1)
                # p = numpy.hstack(results)
                # mean = statistics.mean(p)
                print "update %s - time: %.2f ms" % (str(updateIdx), round(1000.0*(var2-var1),2))
                # print('Mean Elapsed time: ', mean)
                # print >> sys.stderr, '\tRX: "%s"' % data
                # ACK = "ACK"
                # sock.sendto(ACK, address)

            else:
                print >> sys.stderr, 'RX: Done', address
                wfd.close()
                break

    finally:
        # Clean up the connection
        sock.close()
