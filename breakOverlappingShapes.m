function [newshapes,uncontained] = breakOverlappingShapes(shapes)
% should be sorted such that shapes only contain earlier shapes
% once a shape is contained, it cannot be contained again

for i = 1:length(shapes)
	[shapes{i}(1,:),shapes{i}(2,:)] = poly2cw(shapes{i}(1,:),shapes{i}(2,:));
end
uncontained = [];
for i = 1:length(shapes)
	contained = []; % keep track of the uncontained shapes contained by this one
	newshapes{i} = shapes{i};
	for j = 1:length(uncontained)
% 		plot(shapes{i}(1,:),shapes{i}(2,:),...
% 			shapes{uncontained(j)}(1,:),shapes{uncontained(j)}(2,:))
		if inpolygon(shapes{uncontained(j)}(1,1),shapes{uncontained(j)}(2,1),...
				newshapes{i}(1,:),newshapes{i}(2,:))
			newshapes{i} = [[newshapes{i}(1,:),fliplr(shapes{uncontained(j)}(1,:)),shapes{uncontained(j)}(1,end),newshapes{i}(1,end)];...
				[newshapes{i}(2,:),fliplr(shapes{uncontained(j)}(2,:)),shapes{uncontained(j)}(2,end),newshapes{i}(2,end)]];
% 			fill(newshapes{i}(1,:),newshapes{i}(2,:),'r')
			contained = cat(1,contained,j);
		end
	end
	uncontained(contained) = [];
	uncontained = cat(2,uncontained,i);
end