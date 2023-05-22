import cv2
import numpy as np
import os



diffL = []

def detect_color_start(image, column, threshold, numPic):
    img = cv2.imread(image)
    img_height, img_width, _ = img.shape

    pixelOI = img[img_height - 100, 150 : 220]
    avg_green = np.mean(pixelOI[:,1])


    # Iterate from the bottom of the column upwards
    for y in range(img_height - 100, 0, -1):
        if avg_green >= threshold:
            # Calculate the position of the marker
            s = (y - 106) * (2 / (1310-106)) # in mm
            # Draw a vertical marker on the image at the detected position
            cv2.line(img, (0, y), (img_width, y), (0, 255, 0), 2)
            cv2.putText(img, f"l={round(s, 4)}e-6 m", (10, y - 10), cv2.FONT_HERSHEY_SIMPLEX, 0.65, (255, 255, 255), 2)
            # Output the detected position
            print(f"Image: {image}, Height: {s}")
            diffL.append(s)
            break
        pixelOI = img[y, 150 : 220]
        avg_green = np.mean(pixelOI[:,1])  # Update average green value




    # Display the image with the marker
    #cv2.imshow("Image with Marker", img)
    cv2.imwrite("Marked_Image_{}.jpg".format(numPic), img)
    cv2.waitKey(0)
    cv2.destroyAllWindows()

def process_images(directory, column, threshold):
    numPic = 0
    # Iterate over all .jpg images in the directory
    for file in os.listdir(directory):
        if file.endswith(".jpg"):
            image_path = os.path.join(directory, file)
            detect_color_start(image_path, column, threshold, numPic)
            numPic += 1

# Specify the directory containing the images
image_directory = "Microfluidics/Pics"

# Specify the column index (0-based) and threshold for green value
column_index = 160
green_threshold = 200

# Process the images in the directory
process_images(image_directory, column_index, green_threshold)
print(diffL)