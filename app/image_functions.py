from PIL import Image, ImageDraw
from shapely.geometry import Polygon, MultiPolygon
import numpy as np
# import cv2
# from skimage.metrics import structural_similarity as ssim
# 
from app.generic_functions import flatten_list#, roof
 
def restrict_image(image: Image.Image, polygon: Polygon | MultiPolygon) -> Image.Image:
  """
  Restricts an image to a polygon.
  ---
  Input:
    - image: An Image (coming from Image.open(BytesIO)) image.
    - polygon: a shapely polygon or multipolygon
  Output:
    - An image of the same type as input.
  """
  x_size, y_size = image.size
  min_lon, min_lat, max_lon, max_lat = polygon.bounds

  mask = Image.new('L', image.size, 0)
  draw = ImageDraw.Draw(mask)
  
  if isinstance(polygon, Polygon):
    polygon_coords = list(polygon.exterior.coords)
  else:
    polygon_coords = flatten_list([list(p.exterior.coords) for p in polygon.geoms])

  transformed_polygon_coords = [(
    (x_size * (lon - min_lon)) // (max_lon - min_lon),
    (y_size * (max_lat - lat)) // (max_lat - min_lat)
    ) for lon, lat in polygon_coords]

  draw.polygon(transformed_polygon_coords, fill=255)
  
  mask_array = np.array(mask)
  
  if image.mode != 'RGBA':
    image = image.convert('RGBA')
  
  image_array = np.array(image)
  
  alpha = mask_array.astype(np.uint8)
  image_array[:, :, 3] = alpha
  
  result = Image.fromarray(image_array)
  return result
# 
# 
# def raster_image(image: Image.Image, factor: int = 2) -> list[Image.Image]:
#   x_size, y_size = image.size
# 
#   image = image.convert("RGB")
#   np_imag = np.asarray(image)
# 
#   images = []
#   for i in range(factor):
#     for j in range(factor):
#       x_step = i * x_size // factor
#       y_step = j * y_size // factor
# 
#       raster_image = Image.fromarray(
#         np_imag[y_step:y_step + roof(y_size / factor), x_step:x_step + roof(x_size / factor)]
#       )
#       images.append(raster_image)
#   return images
# 
# 
# """ COMPLETELY BORKED! """
# def fill_crop_to_size(image: Image.Image | np.ndarray, target_size: tuple[int, int]) -> Image.Image | np.ndarray:
#     """
#     Resize image to target size while maintaining aspect ratio and adding black padding if necessary.
#     
#     Args:
#         image: PIL.Image.Image or numpy.ndarray
#         target_size: tuple of (width, height)
#         
#     Returns:
#         Same type as input (PIL.Image.Image or numpy.ndarray)
#     """
#     # Convert numpy array to PIL Image if necessary
#     input_is_numpy = isinstance(image, np.ndarray)
#     if input_is_numpy:
#       if image.shape[0] == target_size[0] and image.shape[1] == target_size[1]:
#         return image
#     else:
#       if image.height == target_size[0] and image.width == target_size[1]:
#         return image
# 
#     if input_is_numpy:
#         if image.dtype != np.uint8:
#             image = (image * 255).astype(np.uint8)
#         pil_image = Image.fromarray(image)
#     else:
#         pil_image = image
# 
#     # Get current dimensions
#     current_width, current_height = pil_image.size
#     target_height, target_width = target_size
#     
#     # Calculate scaling factors
#     width_ratio = target_width / current_width
#     height_ratio = target_height / current_height
#     
#     # Use the smaller ratio to ensure the image fits within target size
#     scale_factor = min(width_ratio, height_ratio)
#     
#     # Calculate new dimensions
#     new_width = int(current_width * scale_factor)
#     new_height = int(current_height * scale_factor)
#     
#     # Resize image
#     resized_image = pil_image.resize((new_width, new_height), Image.Resampling.LANCZOS)
#     
#     # Create new black image with target dimensions
#     result_image = Image.new('RGB', (target_width, target_height), 'black')
#     
#     # Calculate position to paste resized image
#     
#     # Paste resized image onto black background
#     result_image.paste(resized_image, (0,0))
#     
#     # Convert back to numpy array if input was numpy
#     if input_is_numpy:
#         return np.array(result_image)
#     return result_image
# 
# def images_similar(img1, img2, mean_threshold=10, sscore_threshold=0.5):
#   """ img1 >= img2 """
# 
#   img1 = np.array(img1.convert("RGB"), dtype=np.uint8)
#   img2 = np.array(img2.convert("RGB"), dtype=np.uint8)
# 
#   # Resize img1 to match img2
#   img1_resized = cv2.resize(img1, (img2.shape[1], img2.shape[0])).astype(np.uint8)
# 
#   # Compute absolute difference
#   diff = np.abs(img1_resized.astype(np.int16) - img2.astype(np.int16))
#   mean_diff = np.mean(diff)
# 
#   # Compute SSIM for RGB
#   sscore = ssim(img1_resized, img2, data_range=255, channel_axis=2)
# 
#   return (mean_diff < mean_threshold) or (sscore > sscore_threshold)
