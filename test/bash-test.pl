




awk '{print "PrediXcan.py --predict --weights "$2" --dosages genotypes --samples genotypes/samples.txt --output_dir outDir/intermediate/"$1" --linear"}' index.txt | xargs -I {} -P $parallel sh -c \"