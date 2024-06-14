import os
from PIL import Image, ImageOps
import math

def get_png_files(directory):
    return [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.png')]

def select_equally_spaced_frames(frames, count):
    step = len(frames) / count
    return [frames[math.floor(i * step)] for i in range(count)]

def add_border(image, border_size=10, color='black'):
    return ImageOps.expand(image, border=border_size, fill=color)

def concatenate_frames_horizontally(frames, border_size=10):
    bordered_frames = [add_border(Image.open(frame), border_size=border_size) for frame in frames]
    widths, heights = zip(*(frame.size for frame in bordered_frames))
    total_width = sum(widths)
    max_height = max(heights)
    
    concatenated_image = Image.new('RGB', (total_width, max_height))
    
    x_offset = 0
    for frame in bordered_frames:
        concatenated_image.paste(frame, (x_offset, 0))
        x_offset += frame.width
    
    return concatenated_image

def process_frames(directory, output1, output2, frame_count=6):
    frames = get_png_files(directory)
    frames.sort()  # Ensure the frames are in order
    
    half_point = len(frames) // 2
    first_half = frames[:half_point]
    second_half = frames[half_point:]
    
    selected_frames_first_half = select_equally_spaced_frames(first_half, frame_count)
    selected_frames_second_half = select_equally_spaced_frames(second_half, frame_count)
    
    concatenated_first_half = concatenate_frames_horizontally(selected_frames_first_half)
    concatenated_second_half = concatenate_frames_horizontally(selected_frames_second_half)
    
    concatenated_first_half.save(output1)
    concatenated_second_half.save(output2)

if __name__ == "__main__":
    directory = './build/frames'
    output1 = './output_first_half.png'
    output2 = './output_second_half.png'
    
    process_frames(directory, output1, output2)
