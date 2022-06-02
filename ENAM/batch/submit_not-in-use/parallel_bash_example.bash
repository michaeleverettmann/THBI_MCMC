#/bin/bash 

task(){
   sleep 0.5; echo "$1";
}

for thing in a b c d e f g; do 
  task "$thing" &
done 

wait
echo "All done"