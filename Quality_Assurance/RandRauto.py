import SimpleITK as sitk
import numpy as np
import tempfile
import nrrd
import os

def Register(pathFixed, pathMoving, pathSegmentation, flagRadiomics=0):

    fixed_image = sitk.ReadImage(pathFixed, sitk.sitkFloat32)
    moving_image = sitk.ReadImage(pathMoving, sitk.sitkFloat32)
    segmentation = sitk.ReadImage(pathSegmentation, sitk.sitkFloat32)

    # Resampling reference image to fixed image
    identity = sitk.Transform(3, sitk.sitkIdentity)
    moving_resampled = sitk.Resample(moving_image, fixed_image, identity, sitk.sitkNearestNeighbor, 0.0, moving_image.GetPixelID())
    segmentation_resampled = sitk.Resample(segmentation, fixed_image, identity, sitk.sitkNearestNeighbor, 0.0, moving_image.GetPixelID())

    initial_transform = sitk.CenteredTransformInitializer(fixed_image, 
                                                      moving_resampled, 
                                                      sitk.Euler3DTransform(), 
                                                      sitk.CenteredTransformInitializerFilter.MOMENTS)

    registration_method = sitk.ImageRegistrationMethod()

    # Similarity metric settings.
    registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    registration_method.SetMetricSamplingStrategy(registration_method.RANDOM)
    registration_method.SetMetricSamplingPercentage(0.01)

    registration_method.SetInterpolator(sitk.sitkLinear)

    # Optimizer settings.
    registration_method.SetOptimizerAsGradientDescent(learningRate=1.0, numberOfIterations=100, convergenceMinimumValue=1e-6, convergenceWindowSize=10)
    registration_method.SetOptimizerScalesFromPhysicalShift()

    # Don't optimize in-place, we would possibly like to run this cell multiple times.
    registration_method.SetInitialTransform(initial_transform, inPlace=False)

    final_transform = registration_method.Execute(fixed_image, moving_resampled)

    # Always check the reason optimization terminated.
    print('Final metric value: {0}'.format(registration_method.GetMetricValue()))
    print('Optimizer\'s stopping condition, {0}'.format(registration_method.GetOptimizerStopConditionDescription()))

    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(fixed_image)
    resampler.SetInterpolator(sitk.sitkLinear)
    resampler.SetDefaultPixelValue(100)
    resampler.SetTransform(final_transform)

    out = resampler.Execute(segmentation_resampled)

    # Transpose 
    mask_test = sitk.GetArrayFromImage(out)
    mask_test[mask_test == 100] = 0
    mask_test = (np.transpose(mask_test, (2,1,0)))

    if flagRadiomics == 1:
        fixedCT, metaFixed = nrrd.read(pathFixed)
        temp_dir = tempfile.TemporaryDirectory()
        nrrd.write(os.path.join(temp_dir.name, 'mask.nrrd'), mask_test, metaFixed)
        print(os.listdir(temp_dir.name))

        return mask_test, temp_dir

    return mask_test

def RegisterOld(pathFixed, pathMoving, pathSegmentation):

    fixed_image = sitk.ReadImage(pathFixed, sitk.sitkFloat32)
    moving_image = sitk.ReadImage(pathMoving, sitk.sitkFloat32)
    segmentation = sitk.ReadImage(pathSegmentation, sitk.sitkFloat32)

    # Resampling reference image to fixed image
    identity = sitk.Transform(3, sitk.sitkIdentity)
    moving_resampled = sitk.Resample(moving_image, fixed_image, identity, sitk.sitkNearestNeighbor, 0.0, moving_image.GetPixelID())
    segmentation_resampled = sitk.Resample(segmentation, fixed_image, identity, sitk.sitkNearestNeighbor, 0.0, moving_image.GetPixelID())

    initial_transform = sitk.CenteredTransformInitializer(fixed_image, 
                                                      moving_resampled, 
                                                      sitk.Euler3DTransform(), 
                                                      sitk.CenteredTransformInitializerFilter.MOMENTS)

    registration_method = sitk.ImageRegistrationMethod()

    # Similarity metric settings.
    registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    registration_method.SetMetricSamplingStrategy(registration_method.RANDOM)
    registration_method.SetMetricSamplingPercentage(0.01)

    registration_method.SetInterpolator(sitk.sitkLinear)

    # Optimizer settings.
    registration_method.SetOptimizerAsGradientDescent(learningRate=1.0, numberOfIterations=100, convergenceMinimumValue=1e-6, convergenceWindowSize=10)
    registration_method.SetOptimizerScalesFromPhysicalShift()

    # Don't optimize in-place, we would possibly like to run this cell multiple times.
    registration_method.SetInitialTransform(initial_transform, inPlace=False)

    final_transform = registration_method.Execute(fixed_image, moving_resampled)

    # Always check the reason optimization terminated.
    print('Final metric value: {0}'.format(registration_method.GetMetricValue()))
    print('Optimizer\'s stopping condition, {0}'.format(registration_method.GetOptimizerStopConditionDescription()))

    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(fixed_image)
    resampler.SetInterpolator(sitk.sitkLinear)
    resampler.SetDefaultPixelValue(100)
    resampler.SetTransform(final_transform)

    out = resampler.Execute(segmentation_resampled)

    # Transpose 
    mask_test = sitk.GetArrayFromImage(out)
    mask_test[mask_test == 100] = 0
    mask_test = (np.transpose(mask_test, (2,1,0)))

    return mask_test

