# A code which downloads the SBND Moliere radius code from the gpvms and adds it to my github page 
scp kmistry@sbndgpvm01.FNAL.GOV:/sbnd/app/users/kmistry/sbndcode_v06_79_00/srcs/sbndcode/sbndcode/MolRadReco/*  .;

git add *;

echo "Now run: git commit -m Enter your update here to save the code" 
echo "Then Run: git push origin master"

