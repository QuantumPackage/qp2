#!/bin/sh   

list='ls *.org'
for element in $list    
do   
		emacs --batch $element -f org-babel-tangle
done
