#! /bin/bash
mkdir -p ./outputs
n=21

for ((i=1;i<=$n;i++)); 
do 
   # your-unix-command-here
    echo $i
    docker stop dispredict_$i && docker rm dispredict_$i
done
sleep 5
echo "All containers are stopped"

for ((i=1;i<=$n;i++)); 
do 
   # your-unix-command-here
    echo $i
    docker stop dispredict_$i && docker rm dispredict_$i
    docker run -itd --name dispredict_$i wasicse/dispredict3.0:latest
    docker cp ./inputs/processedinput_$i.fasta dispredict_$i:/opt/Dispredict3.0/example/sample.fasta
    sleep 2
    docker exec -it dispredict_$i /bin/bash -c "source /opt/Dispredict3.0/.venv/bin/activate && /opt/Dispredict3.0/.venv/bin/python /opt/Dispredict3.0/script/Dispredict3.0.py -f /opt/Dispredict3.0/example/sample.fasta -o /opt/Dispredict3.0/output/" | tee -a ./outputs/output_$i.txt &

done
echo "Waiting for the output"
wait


echo "Output is ready"
for ((i=1;i<=$n;i++)); 
do 
    mkdir -p ./outputs/output_$i
    docker cp dispredict_$i:/opt/Dispredict3.0/output/ ./outputs/output_$i
done
echo "Output is copied"