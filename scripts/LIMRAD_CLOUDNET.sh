#!/bin/bash
#python LIMRAD94_to_Cloudnet.py 20181127 000000 240000
#python LIMRAD94_to_Cloudnet.py 20181128 000000 240000
#python LIMRAD94_to_Cloudnet.py 20181129 000000 240000
#python LIMRAD94_to_Cloudnet.py 20181130 000000 240000

# november
#for i in {27..30}
#do
#	python LIMRAD94_to_Cloudnet.py 201811$i 000000 240000
#done

# dezember
#for i in {01..31}
#do
#	python LIMRAD94_to_Cloudnet.py 201812$i 000000 240000
#done

# january
for i in {21..31}
do
	python LIMRAD94_to_Cloudnet.py 201901$i 000000 240000
done

# february
#for i in {01..03}
#do
#	python LIMRAD94_to_Cloudnet.py 201902$i 000000 240000
#done

