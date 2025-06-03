sudo sysctl -w net.inet.tcp.delayed_ack=0
sudo sysctl -w net.inet.tcp.sendspace=128000
sudo sysctl -w net.inet.tcp.recvspace=128000
cd ~/Desktop/java-source/vision4
sudo nice -n -15 java -Xms50m -Xmx1024m -classpath classes:lib/colt.jar:lib/freehep-base.jar:lib/freehep-graphics2d.jar:lib/freehep-graphicsio.jar edu.ucsc.neurobiology.vision.Vision