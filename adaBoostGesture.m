% Written by Robert Lott for Rubine Classification Project
classdef  adaBoostGesture < handle
    %adaBoostGesture Will hold and process values directly related to
    %classification of data samples
    %   Ground truth is assembled from a simple assignment of "class" which
    %   acts as an ID for the ground truth, it should be between 0 and 9 at
    %   present
    
    properties
        %For Classifiers
        weights
        error
        importance
        mu
        %For Template
        class %Classification ID related to what it represents
        samples
        featureMean
        sampleEstimateCovariance %variance
        variance %also variance
        featureCount %how many features we're using
        featureSet %the actual features from all the samples
        covariance %mean covariance of all classes
        
        %For Drawn Object
        groundTruth     % The class identity of the sample, optional arg
        ID              % Sample ID of the sample, optional arg
        features        % feature values of the shape
        oldX            % x position of points
        oldY            % y position of points
        X               % interpolated x points
        Y               % interpolated y points
        center          % center of adjusted image
        t               % time of when each point ws placed
        boundingBox     % 2 opposite corners of bounding box
        
        %         For features
        theta           %angle of points
        dx              %difference in position where points were placed
        dy              %difference in position where points were placed
        dt              %difference in time between points placed
        euDist          %euclidian distance along points
        %         For Prediction
        label
        LikelyLabels
        LikelyHypothesis
    end
    
    methods
        
        function adaBoostGesture = adaBoostGesture(startEmpty,points,groundTruth,ID)
            %             and if i did, it'll check to see if i wanted to just make an
            %             empty one
            if(exist('startEmpty','var') && startEmpty == 1)
                if(exist('groundTruth', 'var'))
                    adaBoostGesture.class = groundTruth;
                end
                %                    adaBoostGesture.loadObj(savedObj)
                return
            else
                startEmpty = 0;
            end
            %             Checks to make sure i didn't call an empty shape
            if ~exist('points','var')
                %           If there aren't any points, you did a bad thing
                error('You forgot some points.')
                return
            end
            if exist('groundTruth','var')
                if ~exist('ID','var')
                    %                     If it doesn't exist, just generate a random ID
                    ID = rand()*50204918;
                end
                adaBoostGesture.storeAsTemplate(groundTruth,ID);
                %                 disp(['Storing class ' , num2str(groundTruth)])
            end
            adaBoostGesture.storeValues(points);
            adaBoostGesture.findDeltas();
            adaBoostGesture.calculateFeatures();
            
            
        end
        function storeAsTemplate(adaBoostGesture,groundTruth,ID)
            %             Why is there an ID? just in case.
            adaBoostGesture.ID = ID;
            adaBoostGesture.groundTruth = groundTruth;
        end
        function [points] = removeClosePoints(adaBoostGesture,points)
            % remove points that are too close to one another that
            % could throw off the processing
            %               badPoints = find(adaBoostGesture.euDist(2:end-1) < 0.01)
            euDist = adaBoostGesture.euDist;
            len=length(euDist);
            badCounter = 1;
            badPoints = 1;
            for i=2:len-2
                if euDist(i) < 0.01
                    euDist(i+1)=euDist(i+1)+euDist(i);
                    badPoints(badCounter )=i;
                    badCounter =badCounter +1;
                end
            end
            if isempty(badPoints)
                return
            end
            if(length(badPoints) > length(euDist(:,1))-3)
                badPoints = badPoints(badPoints < length(euDist(:,1))-3);
            end
            euDist(badPoints,:)=[];
            points(badPoints+2,:) = [];
        end
        
        function findBoundingBoxAndCenter(adaBoostGesture)
            
            minX = min(adaBoostGesture.X);
            minY = min(adaBoostGesture.Y);
            maxX = max(adaBoostGesture.X);
            maxY = max(adaBoostGesture.Y);
            adaBoostGesture.boundingBox = [minX,minY; maxX,maxY];
            
        end
        
        function findCenter(adaBoostGesture)
            meanX = mean(adaBoostGesture.X );
            meanY = mean(adaBoostGesture.Y);
            adaBoostGesture.center = [meanX, meanY];
        end
        function calculateDistanceBetweenPoints(adaBoostGesture,points)
            adaBoostGesture.euDist = sqrt(...
                sum((points(3:end,1:2)...
                - points(1:end-2,1:2)).^2,2));
        end
        
        function adaBoostGesture = storeAsSample(adaBoostGesture,groundTruth,ID)
            adaBoostGesture.ID = ID;
            adaBoostGesture.groundTruth = groundTruth;
        end
        function adaBoostGesture = storeValues(adaBoostGesture,points)
            adaBoostGesture.euDist = [0;sqrt(sum((points(2:end,1:2) - points(1:end-1,1:2)).^2,2))];
            points = adaBoostGesture.removeClosePoints(points);
            %             points = adaBoostGesture.smoothData(points,5)
            points(:,1) = adaBoostGesture.smoothData(points(:,1),5);
            points(:,2) = adaBoostGesture.smoothData(points(:,2),5);
            points(:,3) = adaBoostGesture.smoothData(points(:,3),5);
            points = points(4:end-4,:);
            adaBoostGesture.euDist = [0;sqrt(sum((points(2:end,1:2) - points(1:end-1,1:2)).^2,2))];
            adaBoostGesture.oldX = points(:,1);
            adaBoostGesture.X= points(:,1);
            adaBoostGesture.oldY = points(:,2);
            adaBoostGesture.Y = points(:,2);
            adaBoostGesture.t = points(:,3);
        end
        function adaBoostGesture = calculateFeatures(adaBoostGesture)
            xMax = max(adaBoostGesture.X);
            yMax = max(adaBoostGesture.Y);
            xMin = min(adaBoostGesture.X);
            yMin = min(adaBoostGesture.Y);
            if(isempty(adaBoostGesture.dx))
                adaBoostGesture.findDeltas();
                disp('Shape is empty, something is wrong.')
            end
            f1 = acos((adaBoostGesture.dx(1)) / sqrt( ...
                (adaBoostGesture.dx(1)).^2 + ...
                (adaBoostGesture.dy(1)).^2));
            f2 = asin((adaBoostGesture.dy(1)) / sqrt( ...
                (adaBoostGesture.dx(1)).^2 +...
                (adaBoostGesture.dy(1)).^2));
            f3 = sqrt((xMax - xMin).^2 + (yMax - yMin).^2);
            f4 = atan2((yMax - yMin) , (xMax - xMin));
            f5 = sqrt((adaBoostGesture.X(end-1) - adaBoostGesture.X(1)).^2 +...
                (adaBoostGesture.Y(end-1) - adaBoostGesture.Y(1)).^2);
            f6 = acos((adaBoostGesture.X(end-1) - adaBoostGesture.X(1))./f5);
            f7 = asin((adaBoostGesture.Y(end-1) - adaBoostGesture.Y(1))./f5);
            f8 = sum(adaBoostGesture.euDist);
            f9 = sum(adaBoostGesture.theta(1:end-2));
            f10 = sum(abs(adaBoostGesture.theta(1:end-2)));
            f11 = sum(adaBoostGesture.theta(1:end-2).^2);
            f12 = max((adaBoostGesture.dx(1:end-2).^2 + ...
                adaBoostGesture.dy(1:end-2).^2)./adaBoostGesture.dt(1:end-2).^2);
            f13 = adaBoostGesture.t(end-1) - adaBoostGesture.t(1);
            adaBoostGesture.features = [...
                f1;f2;f3;f4;f5;f6;f7;...
                f8;f9;f10;f11;f12;f13];
        end
        function adaBoostGesture = findDeltas(adaBoostGesture)
            % calculate the distance between every two data points
            adaBoostGesture.dx = adaBoostGesture.X(3:end)- adaBoostGesture.X(1:end-2);
            adaBoostGesture.dy = adaBoostGesture.Y(3:end)- adaBoostGesture.Y(1:end-2);
            
            adaBoostGesture.dt = adaBoostGesture.t(3:end)-adaBoostGesture.t(1:end-2);
            adaBoostGesture.theta = ...
                atan2(adaBoostGesture.dx(3:end).*adaBoostGesture.dy(1:end-2)...
                - adaBoostGesture.dx(1:end-2).*adaBoostGesture.dy(3:end),...
                adaBoostGesture.dx(3:end).*adaBoostGesture.dx(1:end-2)...
                +adaBoostGesture.dy(3:end).*adaBoostGesture.dy(1:end-2));
        end
        function weight = polarWeight(adaBoostGesture)
            %         A normalization for r values so further away values don't have
            %         as strong an impact by taking the 10th root
            weight =  adaBoostGesture.r.^(0.10);
        end
        function createThetaVarients(adaBoostGesture)
            len = 20;
            shifts = linspace(-pi,pi,len);
            thetaVar = (ones(length(adaBoostGesture.theta),len) .* adaBoostGesture.theta)+ shifts;
            
            adaBoostGesture.thetaVarients = thetaVar;
        end
        function saveStruct = saveobj(adaBoostGesture)
            %             loads object from save file
            saveStruct.groundTruth = adaBoostGesture.groundTruth;
            saveStruct.ID = adaBoostGesture.ID;
            saveStruct.oldX = adaBoostGesture.oldX;
            saveStruct.oldY = adaBoostGesture.oldY;
            saveStruct.X = adaBoostGesture.X;
            saveStruct.Y = adaBoostGesture.Y;
            saveStruct.t = adaBoostGesture.t;
            saveStruct.center = adaBoostGesture.center;
            saveStruct.theta = adaBoostGesture.theta;
            saveStruct.features = adaBoostGesture.features;
        end
        %%%%%%%%%%%%% CLASS GROUPING FUNCTIONS %%%%%%%%%%%%%%%%
        function adaBoostGesture = makeClassifier(adaBoostGesture, adaBoostGestureectArray)
            adaBoostGesture.storeInput(adaBoostGestureectArray);
            adaBoostGesture.fmean();
            adaBoostGesture.storeFeatureCounts();
            adaBoostGesture.calculateSampleEstimateCovariance();
            
            
        end
        
        function adaBoostGesture = storeInput(adaBoostGesture, samples)
            adaBoostGesture.samples = samples;
            featureSet = samples(1).features;
            for i = 2:length(samples)
                featureSet = [featureSet,samples(i).features];
            end
            adaBoostGesture.featureSet = featureSet;
            
        end
        
        function adaBoostGesture = fmean(adaBoostGesture)
            %             finds the mean
            adaBoostGesture.featureMean = mean(adaBoostGesture.featureSet,2);
            
        end
        function adaBoostGesture = storeFeatureCounts(adaBoostGesture)
            %             Instantiates variables for use and stores number of features
            adaBoostGesture.featureCount = length(adaBoostGesture.featureMean);
            adaBoostGesture.covariance = zeros(adaBoostGesture.featureCount,1);
            adaBoostGesture.sampleEstimateCovariance = ...
                zeros(adaBoostGesture.featureCount,adaBoostGesture.featureCount);
        end
        function adaBoostGesture = calculateSampleEstimateCovariance(adaBoostGesture)
            %             Calculates sample estimate covariance
            adaBoostGesture.sampleEstimateCovariance = ...
                (adaBoostGesture.featureSet - adaBoostGesture.featureMean)...
                *(adaBoostGesture.featureSet - adaBoostGesture.featureMean)';
            adaBoostGesture.variance = ...
                sum(adaBoostGesture.sampleEstimateCovariance,2);
            
        end
        function saveStruct = saveClassifier(adaBoostGesture)
            %             loads object from save file
            adaBoostGestureectSamples = adaBoostGesture.samples;
            for i = 1:length(adaBoostGestureectSamples)
                saveStruct.samples(i) = adaBoostGestureectSamples(i).saveobj();
                
                saveStruct.class = adaBoostGesture.class; %Classification ID related to what it represents
                saveStruct.featureMean = adaBoostGesture.featureMean;
                saveStruct.featureSet = adaBoostGesture.featureSet;
                saveStruct.covariance = adaBoostGesture.covariance;
                saveStruct.featureCount = adaBoostGesture.featureCount;
                saveStruct.weights = adaBoostGesture.weights;
                saveStruct.mu = adaBoostGesture.mu;
                
            end
            
        end
        %%%%%%%%%%%%% END CLASS GROUPING FUNCTIONS %%%%%%%%%%%%%%%%
    end
    
    
    methods(Static)
        function smoothedData = smoothData(data,windowSize)
            smoothedData = data;
            len = length(data);
            % Create the kernel
            %             kernel = smoothKernel();
            kernel = [1;4;6;4;1];
            % pad the data so nothing is lost during smoothing
            padding = ceil(windowSize./2);
            edge = padding-1;
            paddedData = zeros(len+2*edge,1);
            paddedData(1+edge:end-edge) = data;
            len2 = len+padding;
            % Now apply the filter
            for i = 1+edge:len2-edge
                smoothedData(i) = sum(paddedData(i-edge:i+edge) .* kernel);
            end
            smoothedData(1) = [];
            
            %                 function kernel = smoothKernel()
            %                 %    Create a gaussian kernel for smoothing
            %                     kernel  = ones(windowSize,1);
            %                     n = windowSize(1)-1;
            %                     for j = 2:windowSize(1)
            %                         kernel (j) = factorial(n)/(factorial(j)*factorial(abs(n-j)));
            %                     end
            %                     kernel = kernel  ./ sum( kernel);
            %                 end
        end
        function loadedObj = loadobj(saved)
            %             loads object from save file
            if isstruct(saved)
                newObj = adaBoostGesture(1);
                newObj.groundTruth = saved.groundTruth;
                newObj.ID = saved.ID;
                newObj.oldX = saved.oldX;
                newObj.oldY = saved.oldY;
                newObj.X = saved.X;
                newObj.Y = saved.Y;
                newObj.t = saved.t;
                newObj.center = saved.center;
                newObj.theta = saved.theta;
                newObj.features = saved.features;
                loadedObj = newObj;
            else
                loadedObj = saved;
            end
        end
        function loadedObj = loadClassifiers(saved)
            if isstruct(saved)
                newClassifier = adaBoostGesture(1);
                newClassifier.class = saved.class;
                newClassifier.featureMean = saved.featureMean;
                newClassifier.featureSet = saved.featureSet;
                newClassifier.covariance = saved.covariance;
                newClassifier.featureCount = saved.featureCount;
                newClassifier.weights = saved.weights;
                newClassifier.mu= saved.mu;
                adaBoostGestureectSamples = saved.samples;
                for i = 1:length(adaBoostGestureectSamples)
                    samples(i) = adaBoostGesture.loadobj(adaBoostGestureectSamples(i));
                end
                newClassifier.samples = samples;
                loadedObj = newClassifier;
            else
                loadedObj = saved;
            end
        end
        function distance = modifiedHausdorff(obj1, obj2)
            %             Hausdorff distance finds the nearest neighbor from one point
            %             in an array to all the other points in the other
            % With Modified Hausdorff you find the weighted mean of all the minimum
            
            polarWeight = obj1.polarWeight();
            %             thetaMod = linspace(-pi,pi,20);
            thetaLen = size(obj1.thetaVarients,2);
            rLen = length(obj1.r);
            %             rLen2 = length(obj2.r);
            thetaDistances = zeros(thetaLen,1);
            for i = 1:thetaLen
                distHolder = zeros(rLen,1);
                for j = 1:rLen
                    %Find distance between a point and all other points
                    R = (obj1.r(j) - obj2.r);
                    %                     Theta = obj1.thetaVarients(j,i) - obj2.theta;
                    %                     Now find the b weighted minimum of the b distance
                    distHolder(j,1) = polarWeight(j) *...
                        min(adaBoostGesture.euclidean(obj1.r(j), obj2.r,...
                        obj1.thetaVarients(j,i),obj2.theta));
                    %                     min(R);
                    %                         min(adaBoostGesture.bDist(R,Theta))
                end
                thetaDistances(i,1) = mean(distHolder(:,1));
                %                 thetaDistances(i,2) = mean(distHolder(:,2));
            end
            distance = max(thetaDistances);
            % distance = mean(distHolder)
        end
        function dist = euclidean(x1,x2,y1,y2)
            %             Distance between points
            dist = sqrt((x2-x1).^2 + (y2-y1).^2);
        end
        function dist = bDist(distribution1,distribution2)
            %             Bhattacharyya
            %             Calculations the discrete probability distrubtions
            dist = -log(sum(sqrt(distribution1'*distribution2))) ;
        end
        function goodArray = getGoodArray(object)
            %           Returns the remaining good features in the event that I do not
            %           actually want to throw away good data
            goodArray = 1:length(object.weights);
            goodArray(object.featuresToRemove) = [];
        end
        function [bestTemplates, scores] = polarRecognize(drawn, templates)
            drawn.createThetaVarients();
            len = length(templates);
            scores = ones(len,1)*99999999;
            for i = 1:len
                scores1(i) = adaBoostGesture.modifiedHausdorff(drawn,templates(i));
                scores2(i) = adaBoostGesture.modifiedHausdorff(templates(i),drawn);
            end
            scores = max([scores1',scores2'],[],2)
            [scores, indices] = sort(scores);
            %             scores = abs(scores);
            
            figure(1)
            imshow(templates(indices(1)).Img)
            title('Best template img')
            figure(2)
            
            plot(templates(indices(1)).r,templates(indices(1)).theta)
            title('Best template polar')
            figure(3)
            
            imshow(drawn.Img)
            title('drawn img')
            figure(4)
            
            plot(drawn.r, drawn.theta)
            title('drawing polar')
            if scores(1)*1.5 < scores(2)
                scores = scores(1);
                bestTemplates = templates(indices(1)).groundTruth;
                return
            end
            bestTemplates = [templates(indices(1)).groundTruth;...
                templates(indices(2)).groundTruth];
            
        end
        
        function classes = classMeans(classes, trainingLength, sampleLimit)
            % T is number of iterations
            % J is weak learners
            % ht is weak hypothesis
            % C is the weak learner
            % Dt is the training weights
            % alpha is the importance
            % before we start, initialize the weights
            % Just the data is being used for training.
            %             normalizationFactor = 1.3;
            sampleLimit = 5;
            classCount = length(classes);
            %             starting off saying they're all weak
            samplesInClass = length(classes(1).samples);
            sampleCount = classCount*samplesInClass;
            featureCount = length(classes(1).samples(1).features);
            weights1 = ones(featureCount,1)*1/featureCount;
            featureArray = zeros(13,10*10);
            labelsArray = ones(13,10);
            set = ones(13,1)
            for i = 1:classCount
                classes(i).weights = weights1;
                featureArray(1:13,10*(i-1)+1 : 10*i) = classes(i).featureSet;
                labelsArray(:,i) = set*(i-1);
            end
            arr1 = [1:classCount];
            arr2 = ones(classCount, trainingLength);
            arr2 = arr2.*arr1';
            labels = ones(featureCount,1);
            f = waitbar(0,"training...")
            remaining  = trainingLength/trainingLength;
            counter = 1;
            
            for t1 = 1:trainingLength*classCount
                
                if counter >= 10
                    counter = 1;
                    remaining = remaining - 1/trainingLength;
                    %                     pause(0.001)
                    waitbar(1-remaining,f, "training...")
                end
                
                t = counter;
                %                create weak classifier ht from Cj using Xi(j),yi and
%                 c = classes(t);
                %                 for i = 1:length(c.samples)
%                 weights = c.weights;
                dist1 = [1:classCount];
                %                 sampleInd = [1:sampleLimit];
                %                 dist1(t) = [];
                % calculate the mus, now find the smallest differences bewteen these and
                % all classes
%                 labels2 = labels*counter;
                
%                 labels2 = labels2;
%                 featureMean1 = c.featureMean;
                for n = 1:numel(dist1(:))
                     
%                     labels3 = labels * (n-1)
%                     label = labels2;
%                     labels2 = labels;
%                     labels2(:,n) = labels(:,n).*-1;
%                     k = dist1(n);
                    c2 = classes(n);
                    labelSet2 = labelsArray;
                    labelSet2(labelSet2 ~= (n-1)) = -1;
                    labelSet2(labelSet2 == (n-1)) = 1;
                   weights2 = c2.weights;
%                     set1 = c2.featureSet;
%                     featureMean2 = c2.featureMean;
                    fmean = c2.featureMean;
%                     mu = sum(featureArray(:,[1:10:100]) .* weights2,3)./sum(weights2);
                    mu = sum(featureArray(:,[1:10:100]) .* weights2,2)./sum(weights2);
                    hypo = sign(mu - fmean);
%                     hypo2 = hypo;
%                     for j = 1:numel(dist1(:))
%                        fmean2 =  classes(n).featureMean;
%                        hypo2(:,:,j) = mu-fmean2;
%                     end
%             mu2 = sum(features(:,[1:10:100]).* weights2,3)./sum(weights2);

%             hyp1 = adaBoostGesture.manhattanDistance(mu1,fmean1);
%             hyp2 = adaBoostGesture.manhattanDistance(mu2,fmean2);
%             weakH = sign(hyp1 - hyp2);

%                     hypo = adaBoostGesture.getWeakHypothesis(featureArray,featureMean1,...
%                         featureMean2, weights, weights2);
%                     
%                     error = sum((hypo < 1).*weights2,2);
                    error = sum((hypo == labelSet2).*weights2,2);
                    error(error<0.5) = error(error<0.5).*-1;
%                     labels2(error >0.49) = -1;
%                     if error < 0.5
%                         labels2 = labels2.*-1;
% % %                         continue
%                     end
                    error( error > 3) = 3;

                    c2.error = error;%(hypo < 1).*weights2;
                    importance = adaBoostGesture.importanceMeasure(error, length(hypo));
                    c2.importance = importance;
                    newWeights = adaBoostGesture.updateWeights(weights2,importance, hypo, labelSet2);
                    classes(n).weights = newWeights;
%                     labels2 = labels'.*-1;
                end
                counter = counter + 1;
                
            end
            pause(0.1)
            waitbar(1,f, "Training Complete!")
            pause(0.2)
            waitbar(0,f, "Storing values...")
            for i = 1:classCount
                classes(i).mu = classes(i).featureMean .* classes(i).weights;
                waitbar(i*1/classCount,f, "Storing values...")
            end
            pause(0.1)
            waitbar(i*1/classCount,f, "Storing values completed")
            pause(0.2)
            waitbar(0,f, "Getting Hypothesis For Samples")
            
            %             final hypothesis;
            LN = zeros([10,10]);
            HN = LN ;
            for i = 1:classCount
                waitbar(i*1/classCount,f, "Getting Hypothesis For Samples...")
                for j = 1:classCount
                    [H, L] = adaBoostGesture.getHypothesis(classes,classes(i).samples(j));
                    classes(i).samples(j).LikelyLabels = L;
                    classes(i).samples(j).LikelyHypothesis = H;
                    LN(i,j) = L;
                    HN(i,j) = H;
                end
            end
            pause(0.2)
            waitbar(1,f, "Hypothesis Complete, goodbye.")
            close(f)
            LN
            HN
        end
        
        function newWeights=  updateWeights(weights,importance,hypothesis, label)
%             weights = model.weights;
%             importance = model.importance;
            %           model.weights = weights.*exp(-importance*label*hypothesis) / sum(weights);
%             hypothesis = mode(hypothesis);
%             newWeights =  weights.*exp(-importance.*label*hypothesis');
            newWeights = weights.*exp(-importance.*label.*hypothesis) ;
%             newWeights(importance >0) = weights(importance >0);
            newWeights = mean(newWeights./sum(newWeights,1),2);
%             newWeights = mean(newWeights,2);
%             newWeights= newWeights./sum(newWeights);
        end
        
        function [likelihood, labels] = getHypothesis(classifiers, sample)
            labels = [1:10];
            % Prepare all mus
            muSets = zeros([length(classifiers(1).mu),10]);
            weights = zeros([length(classifiers(1).mu),10]);
            fmean = weights;
            for i = 1:length(classifiers)
                classifier = classifiers(i);
                labels(i) = classifier.class;
                muSets(:,i) = classifier.mu;
                fmean(:,i) = classifier.featureMean;
                weights(:,i) = classifier.weights;
            end
            %             construct a distance network for this stuff
            mu = fmean.*weights;
%             distanceSets = adaBoostGesture.manhattanDistance(mu,sample.features);
            distanceSets = mu - sample.features;
            distanceSets2 = distanceSets(:,1);
            indexSets2 = distanceSets2;
            %             indexSets2 = distanceSets;
            for i = 1:13
                [minVal,indices] = min(abs(distanceSets(i,:)));
                indexSets2(i,:) =indices;
                distanceSets2(i,:) = minVal;
                %                 sets2 = sum(sign(distanceSets - distanceSets(:,i)));
                %                distanceSets2(:,i) = sets2;
            end
            %                         distanceSets2 = zeros(10);
            %             distLength = length(distanceSets);
            %             for i = 1:10
            %                 sets2 = sum(sign(distanceSets - distanceSets(:,i)));
            %                distanceSets2(`:,i) = sets2;
            %             end
            %             for i= 1:10
            [labels1, likelihood] = mode(indexSets2);
            likelihood = likelihood/10;
            labels = labels(labels1);
            %                 [labels2,hypothesis2] = mode(indexSets2);
            %                 labels = [labels;labels2];
            %                 hypothesis = [hypothesis;hypothesis2];
            %             end
            %             [hypothesis,index] = sort(sum(distanceSets2));
            %             labels = labels(index);
        end
%         function weakH = getWeakHypothesis(features, fmean1, fmean2,weights1, weights2)
%             mu1 = sum(features .* weights1,2)./sum(weights1);
%             mu2 = sum(features.* weights2,2)./sum(weights2);
% %             mu1 = sum(features .* weights1,3)./sum(weights1);
% %             mu2 = sum(features.* weights2,3)./sum(weights2);
% 
%             hyp1 = adaBoostGesture.manhattanDistance(mu1,fmean1);
%             hyp2 = adaBoostGesture.manhattanDistance(mu2,fmean2);
%             weakH = sign(hyp1 - hyp2);
%         end
         function weakH = getWeakHypothesis(features, fmean, weights, labels)
%             mu1 = sum(features .* weights1,2)./sum(weights1);
%             mu2 = sum(features.* weights2,2)./sum(weights2);
            mu1 = sum(features(:,[1:10:100]) .* weights1,3)./sum(weights1);
%             mu2 = sum(features(:,[1:10:100]).* weights2,3)./sum(weights2);

            hyp1 = adaBoostGesture.manhattanDistance(mu1,fmean1);
%             hyp2 = adaBoostGesture.manhattanDistance(mu2,fmean2);
            weakH = sign(hyp1 - hyp2);
        end
        function dist = euclideanDistance(vector1, vector2)
           dist = sqrt(sum((vector1 - vector2).^2,2));
        end
        function dist = manhattanDistance(vector1,vector2)
            dist = abs(vector1-vector2);
        end
        function  average = weightedAverage(features, weights)
            average = features * weights./ sum(weights);
        end
        function alpha = importanceMeasure(error, labelCount)
            alpha = real((log((1-error) ./ error))); %+ log(labelCount-1);
            alpha(isinf(alpha)) = 0;
        end
        %         function [result,B]=majority_vote(A,weight,numclass)
        %         % FUNCTION: majority voting with weights
        %         % obsv [nr*nc,nbands] each column refers to one feature and each row means one observation
        %         % weight : respective weights of all features (total: nbands)
        %         % numclass: the number of class
        %
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         [d1,d2]=size(A);
        %         B=zeros([d1,numclass]); %weight matrix for each class
        %         for i=1:d2
        %         pos=sub2ind(size(B),1:d1,A(:,i)'); %the value of each column in A corresponds to each class in B
        %         B(pos)=B(pos)+weight(i);
        %         end
        %         [maxv,result]=max(B,[],2); %
        %         total=sum(B==maxv,2);%
        %         %IF there exists multiple maximum, then result will choose the smallest value in A.
        %         unvalid=0;
        %         if sum(total(:)>1)
        %         disp('exist some multiple maximum');
        %         end
        %     end
    end
end

