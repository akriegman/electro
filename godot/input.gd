extends Node3D

signal moved(Vector3)
signal made(Vector3)

@export var joystick: VirtualJoystick
@export var dampening = 0.1
@export var force = 4.0
@export var sensitivity = 200
var velocity = Vector3.ZERO

func _ready():
	Input.set_mouse_mode(Input.MOUSE_MODE_CAPTURED)
	made.emit(position)

func _process(delta):
	var acceleration = Vector3.ZERO
	
	if Input.is_action_pressed("move_left"):
		acceleration += Vector3.LEFT
	if Input.is_action_pressed("move_right"):
		acceleration += Vector3.RIGHT
	if Input.is_action_pressed("move_forward"):
		acceleration += Vector3.FORWARD
	if Input.is_action_pressed("move_back"):
		acceleration += Vector3.BACK
	if Input.is_action_pressed("move_up"):
		acceleration += Vector3.UP
	if Input.is_action_pressed("move_down"):
		acceleration += Vector3.DOWN
	
	if joystick.is_pressed:
		acceleration += \
			joystick.output.x * Vector3.RIGHT + \
			joystick.output.y * Vector3.BACK
	
	acceleration = acceleration.normalized().rotated(Vector3.UP, rotation.y) * force
	velocity += acceleration * delta
	velocity *= dampening ** delta
	position += velocity * delta
	moved.emit(position)
	
func _unhandled_input(event):
	if event is InputEventMouseMotion:
		rotation.y -= event.relative.x / sensitivity
		rotation.x -= event.relative.y / sensitivity

func _on_up_button_down():
	Input.action_press("move_up")

func _on_up_button_up():
	Input.action_release("move_up")

func _on_down_button_down():
	Input.action_press("move_down")

func _on_down_button_up():
	Input.action_release("move_down")
