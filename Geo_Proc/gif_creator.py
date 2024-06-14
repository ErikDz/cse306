from PIL import Image
import os

def create_gif(input_folder, output_gif, duration):
    # Get all frame file names and sort them
    frames = [f for f in os.listdir(input_folder) if f.endswith('.png')]
    frames.sort(key=lambda f: int(f.split('_')[1].split('.')[0]))

    # Load all frames into a list
    images = [Image.open(os.path.join(input_folder, frame)) for frame in frames]

    # Save as GIF
    images[0].save(output_gif, save_all=True, append_images=images[1:], duration=duration, loop=0)

# Usage example
input_folder = 'build/frames'
output_gif = 'n100.gif'
duration = 100  # Duration of each frame in milliseconds

create_gif(input_folder, output_gif, duration)
