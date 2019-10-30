========================
dcicutils for tibanna_ff
========================

Import json to manipulate json format and ff_utils from the dcicutils package.

.. code-block:: python

    import json
    from dcicutils import ff_utils


Key and Authentication
^^^^^^^^^^^^^^^^^^^^^^

S3_ENCRYPT_KEY must be set for the environment.

.. code-block:: bash

    export S3_ENCRYPT_KEY=<s3_encrypt_key>

To obtain the key.

.. code-block:: python

    key = ff_utils.get_authentication_with_server(ff_env='<environment>')

List of environments:

  - fourfront-webdev (4DN portal)
  - fourfront-cgap (CGAP portal)


Posting objects on the portal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Objects are defined with a json structure. Each object type is associated to a specific structure.

1) Reading the object from file.

.. code-block:: python

    fi = open('<path/to/json object>')

2) Loading the object.

.. code-block:: python

    object = json.load(fi)

3) Posting the json object.

.. code-block:: python

    ff_utils.post_metadata(object, '<object type>', key=key)

List of object types:

  - Workflow
  - Software
  - FileReference
  - FileProcessed


Patching objects on the portal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Objects can be patched and updated on the portal. The steps are the same as for posting, but the uuid (unique universal identifier) is required for patching.

.. code-block:: python

    ff_utils.patch_metadata(object_updated, '<uuid>', key=key)


Accessing objects metadata
^^^^^^^^^^^^^^^^^^^^^^^^^^

Metadata can be accessed using the uuid associated with the object. Metadata are returned as dict.

.. code-block:: python

    dict_metadata = ff_utils.get_metadata('<uuid>', key=key)

.. note::

  E.g. to access quality metrics info for a fastq file we can use ff_utils.get_metadata() to fetch the fastq object metadata using its uuid. From the object dict we can then retrieve the uuid for the metrics object and use it to fetch the metrics we need.
