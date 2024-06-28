
function h=mnormalize(u);
up=u-min(u(:));
h=up./max(up(:));
