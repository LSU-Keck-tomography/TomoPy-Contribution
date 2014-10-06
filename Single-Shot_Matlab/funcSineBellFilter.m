function [image]=funcSineBellFilter(image)
%[rows,columns,deltaRow,sineBellRow,deltaColumn,sineBellColumn,filteredImage]=funcSineBellFilter(image)
[rows,columns]=size(image);
filteredImage=zeros([rows,columns]);
deltaRow=1/(rows+1);
sineBellRow=sin(pi*(deltaRow:deltaRow:(1-deltaRow)));
deltaColumn=1/(columns+1);
sineBellColumn=sin(pi*(deltaColumn:deltaColumn:(1-deltaColumn)));
for r=1:rows
    for c=1:columns
    filteredImage(r,c)=image(r,c)*sineBellRow(1,r)*sineBellColumn(1,c);
    end
end
image=filteredImage;
