    1  sudo apt update
    2  sudo apt upgrade -y
    3  sudo apt install -y r-base
    4  cat /etc/os-release 
    5  sudo apt install --no-install-recommends software-properties-common dirmngr
    6  wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
    7  sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
    8  sudo apt install --no-install-recommends r-base
    9  sudo apt-get install gdebi-core
   10  wget https://download2.rstudio.org/server/focal/amd64/rstudio-server-2024.12.1-563-amd64.deb
   11  sudo gdebi rstudio-server-2024.12.1-563-amd64.deb
   12  sudo systemctl status rstudio-server
   13  sudo systemctl start rstudio-server
   14  sudo systemctl status rstudio-server
   15  sudo nanp /etc/rstudio/rserver.conf
   16  sudo nano /etc/rstudio/rserver.conf
   17  less /etc/rstudio/rserver.conf
   18  sudo systemctl restart rstudio-server
   19  curl http://localhost:8787
   20  sudo journalctl -u rstudio-server
   21  sudo passwd lirongp
   22  exit
   23  pwd
   24  df -h
   25  gsutil ls gs://h11-sc-bucket
   26  git clone https://github.com/plrlhb12/sc_mouse_bladder.git
   27  ls -lth
   28  cd sc_mouse_bladder/
   29  ls -lth
   30  ls notebooks/
   31  ls scripts/
   32  mkdir -p data results
   33  ls -lth
   34  less .gitignore 
   35  vim .gitignore 
   36  less .gitignore 
   37  git log
   38  git add .
   39  git commit -m "modify .gitignore to keep empty data and results folder"
   40  git config --global user.email "lrpeng@hotmail.com"
   41  git config --global user.name "plrlhb12"
   42  git commit -m "modify .gitignore to keep empty data and results folder"
   43  git push
   44  ls ~/.ssh/id_rsa.pub
   45  touch ~/.ssh/id_rsa.pub
   46  ssh-keygen -t ed25519 -C "your_email@example.com"
   47  ssh-keygen -t ed25519 -C "lrpeng@hotmail.com"
   48  cat ~/.ssh/id_rsa.pub
   49  cat /home/lirongp/.ssh/id_ed25519.pub
   50  ssh -T git@github.com
   51  git push
   52  git add remote origin 
   53  git remote -v
   54  git push
   55  git remote set-url origin git@github.com:plrlhb12/sc_mouse_bladder.git
   56  git push
   57  gsutil cp gs://your-bucket-name/yourfile.fastq.gz /home/lirongp/
   58  gsutil cp gs://h11-sc-bucket/* data/
   59  gsutil cp -r gs://h11-sc-bucket/* data/
   60  ls data/*
   61  ls -lth
   62  touch notebooks/gcp_vm_setup.md
   63  history
   64  history >> notebooks/gcp_vm_setup.md 
