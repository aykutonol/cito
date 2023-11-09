import numpy as np

def save_scene(X_desired, file_number):
    # X_desired is a 1D numpy array
    # file_number is an integer

    # Save scene to yaml file
    file_name = "scenes/" + str(file_number) + ".yaml"
    with open(file_name, 'w') as f:
        f.write("desiredFinalPos: " + str(X_desired) + "\n")
        

if __name__ == "__main__":
    num_scenes = 100
    save_location = "data/scenes/"

    desiredFinalPos: [1.3, 0, 0.0646922, 0.0, 0.0, 0.0]

    lower_x = 1.2
    upper_x = 1.4
    lower_y = -0.1
    upper_y = 0.1

    for i in range(num_scenes):
        x = np.random.uniform(lower_x, upper_x)
        y = np.random.uniform(lower_y, upper_y)

        desiredFinalPos = [x, y, 0.0646922, 0.0, 0.0, 0.0]

        save_scene(desiredFinalPos, i)