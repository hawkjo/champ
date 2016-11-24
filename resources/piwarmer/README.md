# Piwarmer Programs

We used [Piwarmers](https://github.com/jimrybarski/piwarmer) for temperature control.

Here, you can find the heating programs we use to perform experiments and regenerate MiSeq chips, as well as the driver - the set of PID parameters that keep the heating blocks at the desired temperature. Precise control of the temperature is critical for the success of CHAMP experiments. 

The programs and driver are respresented as serialized Django models, which is what Piwarmers use internally. To import and save the objects, you can run the code below. This is mostly supplied for purposes of reproducibility and transparency - in reality, it would probably be faster and easier to just use the Piwarmer GUI to create them manually.

WARNING: The PID values in the driver are only valid for the specific heating blocks that we use. Blocks with any other shape, size or material may become extremely hot, up to several hundred degrees Celcius! 

You will need to install some dependencies first if you don't already have them:

```
sudo apt-get install libyaml-dev
sudo pip install pyyaml

```
Then run this on the Piwarmer in the `resources` directory, after cloning this repo:

```
from django.core import serializers

for filename in ("programs.yml", "driver.yml"):
    with open(filename) as f:
        for obj in serializers.deserialize("yaml", f):
            obj.save()

```
