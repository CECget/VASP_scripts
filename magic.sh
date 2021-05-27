#!/bin/bash

####### author: QYY @XJTU ######
####### e-mail: cba19971011@qq.com ####
####### reference:https://www.jb51.net/article/48832.htm ####


# function read_dir(){
# for file in `ls $1` #read the outside input $1
# do
#  if [ -d $1"/"$file ] #注意此处之间一定要加上空格，否则会报错
#  then
#  read_dir $1"/"$file
#  else
#  echo $1"/"$file #在此处处理文件即可
#  fi
# done
# } 
# #读取第一个参数
# read_dir $1

function read_dir(){
ncount=0
for subfile in `ls $1` #read the outside input $1
do
 if [ -d $1"/"$subfile ] #注意此处之间一定要加上空格，否则会报错
 then
 read_dir $1"/"$subfile
 else
 ncount+=1
  if (( $ncount == 1 )) 
  then
  echo $1
  cd $1"/" #此处处理文件
  cifout.sh #此处处理文件
  cd - >/dev/null
  fi
 fi
done
} 
#读取第一个参数
read_dir $1
