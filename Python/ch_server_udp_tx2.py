import socket
import sys
import time
import struct
import statistics
import math

nTxAntennas = 1

# Create UDP socket 
sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
server_address = ('', 50002)
sock.bind(server_address)
print("starting up on IP ", str(server_address[0]), " and port ", str(server_address[1]))

# Multicast communication
# multicast_group = '224.3.29.71'
# group = socket.inet_aton(multicast_group)
# mreq = struct.pack('4sL', group, socket.INADDR_ANY)
# sock.setsockopt(socket.IPPROTO_IP, socket.IP_ADD_MEMBERSHIP, mreq)

buffMax = 8*(2*nTxAntennas + 1)  # Define packet segmentation (characters in TX/RX strings)
name = "weights_tx2.bin"

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
                wfd.seek(0)
                wfd.write(data)
                print("Received data: ", data.encode('hex_codec'))
                var2 = time.time()
                updateIdx = updateIdx + 1;
                print("update ", str(updateIdx) , "- time: ", round(1000.0*(var2-var1),2))
            else:
                print("RX: Done ", address[1])
                wfd.close()
                break
    finally:
        # Clean up the connection
        sock.close()
