function plot_svpos(svpos, userpos, prnstr)

% INPUTS:
%   svpos{NDATAxNSV} (1x3)
%   userpos, meters, ecef
%   prnstr: cell of PRN strings

sz = size(svpos);
nsv = sz(2);
ndt = sz(1);

% restructure
%   svpos (3xNSVxNDATA) , meters, ecef
svpos_ = zeros(3,nsv,ndt);
for s = 1:nsv
  for k = 1:ndt
    svpos_(:,s,k) = svpos{k,s};
  end
end


r_earth = 6371e3;
figure;
[xe,ye,ze] = sphere(20);
surf(...
  xe*r_earth,...
  ye*r_earth,...
  ze*r_earth...
)
xlabel('X'); ylabel('Y'); zlabel('Z');
title('SV Positions over time')
hold on

ax_sz =  max(max(max(abs(svpos_))));
axis equal
xlim([-ax_sz ax_sz]);
ylim([-ax_sz ax_sz]);
zlim([-ax_sz ax_sz]);

plot3(userpos(1),userpos(2),userpos(3), 'k.', 'MarkerSize', 20)

colors = distinguishable_colors(nsv);


for k = 1:nsv
  sat_hist = reshape(svpos_(:,k,:),3,ndt);
  plot3(...
    sat_hist(1,:),sat_hist(2,:),sat_hist(3,:),...
    '.','MarkerSize', 14, 'MarkerEdgeColor', colors(k,:) ...
  )
end



leg = {'Earth','User',prnstr{:}};
legend(leg)

end