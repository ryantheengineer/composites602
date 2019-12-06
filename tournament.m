function [winner] = tournament(parents,size)
% size = number of parents to compete with each other
numParents = length(parents);
for j = 1:numParents
    for i = 1:size
        designNum(i) = 1 + floor(rand*numParents);
        Score(i) = maximin(designNum(i),parents);
    end
    [~, index] = min(Score);
    winner(j) = parents(designNum(index));
end

