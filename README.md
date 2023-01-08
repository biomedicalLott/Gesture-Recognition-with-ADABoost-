# Gesture-Recognition-with-ADABoost-
<p><img src = "https://i.imgur.com/9DrXPNa.png"/> </p>

An implementation of ADABoost in matlab to recognize hand written gestures. The ultimate goal being to implement a robust ensemble method of gesture recognition. 

As per the paper[^1] I would implement a handwriting symbol classifier using features taken from Rubine[^2] and Speed Seg[^3]. 43 features in total, this model takes a long time to train but once trained it moved fairly quickly. Although in its current implementation the ADABoost underperforms I believe with some minor modifications I could have it running at a much better pace. 

The app is very simple, and gestures area already pre-trained.  Simply draw your digit and select "classify". 


*note:* this is very likely broken as it is an earlier version of the one I submitted as I have misplaced the working version. Use at your own risk


[^1]: LaViola, Joseph J., and Robert C. Zeleznik. "A practical approach for writer-dependent symbol recognition using a writer-independent symbol recognizer." IEEE Transactions on pattern analysis and machine intelligence 29.11 (2007): 1917-1926.
[^2]: Rubine, Dean Harris. The automatic recognition of gestures. Carnegie Mellon University, 1991.
[^3]: Herold, James, and Thomas F. Stahovich. "Speedseg: A technique for segmenting pen strokes using pen speed." Computers & Graphics 35.2 (2011): 250-264.


