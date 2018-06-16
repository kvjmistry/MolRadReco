# A code which downloads the LArIAT Moliere radius code from the gpvms and adds it to 
scp kmistry@uboonegpvm01.FNAL.GOV:/uboone/app/users/kmistry/LArIAT_MolRadReco/srcs/lariatsoft/LArIATAnaModule/MolRadReco_module.cc  .;
scp kmistry@uboonegpvm01.FNAL.GOV:/uboone/app/users/kmistry/LArIAT_MolRadReco/srcs/lariatsoft/LArIATAnaModule/molradrecoorig.fcl  .;

git add *;

echo "Now run: git commit -m Enter your update here to save the code" 
echo "Then Run: git push origin master"

