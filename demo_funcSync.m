%Load the data
[Data] = load_dataset('data\Fourleg\');
%Load parameters
load('data\Parameters.mat');
%Compute basis
for i = 1:length(Data.shapes)
  Data.basis{i} = cotangent_basis(Data.shapes{i}, 32);
end
% Optimize consistent functional maps
[Data.opt_fmaps] = joint_fmap_opt_lb(Data, Para);
% Convert functional maps into point-to-point maps
[Data.opt_maps] = batch_func_2_point(Data, Data.opt_fmaps);
% Evaluate the correspondences
[curve] = eval_point_maps(Data, Data.opt_maps, 256);
plot(curve(1:64), 'r');
hold on;
[curve] = eval_point_maps(Data, Data.initial_maps, 256);
plot(curve(1:64), 'b');