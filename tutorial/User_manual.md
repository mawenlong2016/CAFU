### CAFU Docker image installation
- Step 1: [Docker installation](https://github.com/cma2015/CAFU/blob/master/tutorial/Docker_installation.md)
- Step 2: CAFU installation from Docker Hub
  ```bash
  # Pull latest version of CAFU from Docker Hub
  $ docker pull malab/cafu
  ```
- Step 3: Qucikly start
  ```bash
  $ docker run -it -p 80:80 malab/cafu bash
  $ cd /home/galaxy
  $ bash run.sh
  ```
  Then you can access CAFU instance via http://localhost:80


### Upload data
#### Download CAFU testing datasets
- Download test data from [CAFU GitHub project](https://github.com/cma2015/CAFU). Click **Clone or download**, and download the ZIP compressed files into you local device and then uncompress it. 

  ![Download data](https://github.com/cma2015/CAFU/blob/master/CAFU_images/1.png)

- For users who installed [Git](https://gist.github.com/derhuerst/1b15ff4652a867391f03), one can use following command to download CAFU project.
  ```bash
  git clone https://github.com/cma2015/CAFU.git
  ```
#### Upload regular file
- Click **Get Data** to upload files (see figure blow).

  ![Upload files](https://github.com/cma2015/CAFU/blob/master/CAFU_images/2.png)

  And then you will see the following interface:

  ![](https://github.com/cma2015/CAFU/blob/master/CAFU_images/3.png)

  Clicking **Choose local file** and select a file you would like to upload (e.g. upload the file in the directory ```./CAFU/test_data/SE RNA-Seq/mapInfoSE```), you will get this interface:
  
  ![Upload regular file](https://github.com/cma2015/CAFU/blob/master/CAFU_images/4.png)
  
  Then click **Start** to start to upload files.

#### Upload collection file
- Click **Get Data** to upload files (see figure blow).

  ![Upload files](https://github.com/cma2015/CAFU/blob/master/CAFU_images/2.png)
  
  And then you will see the following interface:
  
  ![Upload files](https://github.com/cma2015/CAFU/blob/master/CAFU_images/5.png)
  
  Seleting a list files to upload as a collection:
  
  ![Upload files](https://github.com/cma2015/CAFU/blob/master/CAFU_images/6.png)

  Click **Start** to upload, after finishing uploading, click **Build** (see figure below):
  
  ![Upload files](https://github.com/cma2015/CAFU/blob/master/CAFU_images/7.png)
  
  Then enter a name for you collection and click **Create list**
  
  ![Upload files](https://github.com/cma2015/CAFU/blob/master/CAFU_images/8.png)
