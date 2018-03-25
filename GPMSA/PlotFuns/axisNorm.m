function axisNorm(mode)
% Tool to set plot axes to the same values.
% axisNorm('start') resets the axis set, starts a new normalization
% axisNorm('add') adds the current axes to the axis set
% axisNorm('set') sets all the axes to the same value,
%                 the value is the outer hull of all axes.

persistent aRec aList

switch mode
case 'start'
  aRec=[];
  aList=[];
case 'add'
  a=axis;
  if isempty(aRec)
    switch(length(a))
    case 4
      aRec=[Inf -Inf Inf -Inf];
    case 6
      aRec=[Inf -Inf Inf -Inf Inf -Inf];
    otherwise
      error('Crazy, crazy axes found by axisNorm');
    end
  end
  aRec(1)=min(aRec(1),a(1));
  aRec(2)=max(aRec(2),a(2));
  aRec(3)=min(aRec(3),a(3));
  aRec(4)=max(aRec(4),a(4));
  if length(aRec)==6
    aRec(5)=min(aRec(5),a(5));
    aRec(6)=max(aRec(6),a(6));
  end
  aList=[aList gca];
case 'set'
  for ii=1:length(aList)
    axis(aList(ii),aRec);
  end
end
